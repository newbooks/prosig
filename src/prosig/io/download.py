from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from time import monotonic
from urllib.error import URLError
from urllib.request import Request, urlopen

MIN_THREADED_BYTES = 50 * 1024 * 1024
DEFAULT_THREADS = 16


class DownloadError(RuntimeError):
    """Raised when a file cannot be downloaded successfully."""


@dataclass(frozen=True)
class DownloadResult:
    url: str
    destination: Path
    bytes_written: int
    content_length: int | None
    threaded: bool


ProgressCallback = Callable[[int, int | None], None]
Clock = Callable[[], float]


@dataclass(frozen=True)
class DownloadMetadata:
    content_length: int | None
    accepts_ranges: bool


@dataclass(frozen=True)
class SingleThreadDownload:
    bytes_written: int
    content_length: int | None


def download_file(
    url: str,
    destination: str | Path,
    *,
    chunk_size: int = 1024 * 1024,
    progress: ProgressCallback | None = None,
    progress_interval_seconds: float = 30.0,
    threads: int = DEFAULT_THREADS,
    timeout: float = 60.0,
    clock: Clock = monotonic,
) -> DownloadResult:
    """Download a URL to a destination path, always replacing the destination.

    The caller is responsible for deciding whether an existing file should be
    overwritten. Once called, this function performs a fresh download into a
    temporary ``.part`` file and atomically replaces the destination only after
    the download completes.
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be greater than zero")
    if threads <= 0:
        raise ValueError("threads must be greater than zero")
    if progress_interval_seconds <= 0:
        raise ValueError("progress_interval_seconds must be greater than zero")

    destination_path = Path(destination)
    part_path = destination_path.with_name(f"{destination_path.name}.part")
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    part_path.unlink(missing_ok=True)

    bytes_written = 0
    content_length = None
    last_progress_at = clock()

    try:
        try:
            metadata = inspect_download(url, timeout=timeout)
        except DownloadError:
            metadata = DownloadMetadata(content_length=None, accepts_ranges=False)
        content_length = metadata.content_length
        threaded = should_use_threaded_download(metadata, threads=threads)
        if threaded:
            try:
                bytes_written = _download_threaded(
                    url,
                    part_path,
                    content_length=metadata.content_length,
                    chunk_size=chunk_size,
                    threads=threads,
                    progress=_throttled_progress(
                        progress,
                        metadata.content_length,
                        progress_interval_seconds,
                        clock,
                        last_progress_at,
                    ),
                    timeout=timeout,
                )
            except (OSError, URLError, DownloadError):
                part_path.unlink(missing_ok=True)
                threaded = False
                single_result = _download_single_threaded(
                    url,
                    part_path,
                    chunk_size=chunk_size,
                    expected_content_length=metadata.content_length,
                    progress=_throttled_progress(
                        progress,
                        metadata.content_length,
                        progress_interval_seconds,
                        clock,
                        last_progress_at,
                    ),
                    timeout=timeout,
                )
                bytes_written = single_result.bytes_written
                content_length = single_result.content_length
        else:
            single_result = _download_single_threaded(
                url,
                part_path,
                chunk_size=chunk_size,
                expected_content_length=metadata.content_length,
                progress=_throttled_progress(
                    progress,
                    metadata.content_length,
                    progress_interval_seconds,
                    clock,
                    last_progress_at,
                ),
                timeout=timeout,
            )
            bytes_written = single_result.bytes_written
            content_length = single_result.content_length
        if content_length is not None and bytes_written != content_length:
            raise DownloadError(
                f"downloaded {bytes_written} bytes but expected "
                f"{content_length}: {url}"
            )
        part_path.replace(destination_path)
    except (OSError, URLError, DownloadError) as error:
        part_path.unlink(missing_ok=True)
        if isinstance(error, DownloadError):
            raise
        message = f"failed to download {url} to {destination_path}"
        raise DownloadError(message) from error

    return DownloadResult(
        url=url,
        destination=destination_path,
        bytes_written=bytes_written,
        content_length=content_length,
        threaded=threaded,
    )


def inspect_download(url: str, *, timeout: float = 60.0) -> DownloadMetadata:
    request = Request(url, headers={"User-Agent": "prosig"}, method="HEAD")
    try:
        with urlopen(request, timeout=timeout) as response:
            status = getattr(response, "status", None)
            if status is not None and status >= 400:
                raise DownloadError(f"download failed with HTTP status {status}: {url}")
            return _metadata_from_headers(response.headers)
    except (OSError, URLError) as error:
        message = f"failed to inspect download metadata for {url}"
        raise DownloadError(message) from error


def should_use_threaded_download(
    metadata: DownloadMetadata,
    *,
    threads: int = DEFAULT_THREADS,
) -> bool:
    if threads <= 1:
        return False
    if not metadata.accepts_ranges:
        return False
    if metadata.content_length is None:
        return False
    return metadata.content_length > MIN_THREADED_BYTES


def _metadata_from_headers(headers: object) -> DownloadMetadata:
    content_length_header = headers.get("Content-Length")
    content_length = int(content_length_header) if content_length_header else None
    accept_ranges = (headers.get("Accept-Ranges") or "").lower()
    return DownloadMetadata(
        content_length=content_length,
        accepts_ranges=accept_ranges == "bytes",
    )


def _download_single_threaded(
    url: str,
    part_path: Path,
    *,
    chunk_size: int,
    expected_content_length: int | None,
    progress: ProgressCallback | None,
    timeout: float,
) -> SingleThreadDownload:
    request = Request(url, headers={"User-Agent": "prosig"})
    bytes_written = 0
    content_length = expected_content_length
    with urlopen(request, timeout=timeout) as response:
        status = getattr(response, "status", None)
        if status is not None and status >= 400:
            raise DownloadError(f"download failed with HTTP status {status}: {url}")
        if content_length is None:
            content_length = _metadata_from_headers(response.headers).content_length

        with part_path.open("wb") as output:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                output.write(chunk)
                bytes_written += len(chunk)
                if progress is not None:
                    progress(bytes_written, content_length)
    return SingleThreadDownload(
        bytes_written=bytes_written,
        content_length=content_length,
    )


def _download_threaded(
    url: str,
    part_path: Path,
    *,
    content_length: int | None,
    chunk_size: int,
    threads: int,
    progress: ProgressCallback | None,
    timeout: float,
) -> int:
    if content_length is None:
        raise DownloadError("threaded download requires a known content length")

    ranges = _byte_ranges(content_length, threads)
    with part_path.open("wb") as output:
        output.truncate(content_length)

    progress_increment = _progress_accumulator(progress)
    with ThreadPoolExecutor(max_workers=len(ranges)) as executor:
        futures = [
            executor.submit(
                _download_range,
                url,
                part_path,
                start,
                end,
                chunk_size,
                timeout,
                progress_increment,
            )
            for start, end in ranges
        ]
        bytes_written = sum(future.result() for future in futures)

    return bytes_written


def _download_range(
    url: str,
    part_path: Path,
    start: int,
    end: int,
    chunk_size: int,
    timeout: float,
    progress_increment: Callable[[int], None] | None,
) -> int:
    request = Request(
        url,
        headers={
            "User-Agent": "prosig",
            "Range": f"bytes={start}-{end}",
        },
    )
    bytes_written = 0
    offset = start
    with urlopen(request, timeout=timeout) as response:
        status = getattr(response, "status", None)
        if status not in (None, 206):
            message = f"range download failed with HTTP status {status}: {url}"
            raise DownloadError(message)
        with part_path.open("r+b") as output:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                output.seek(offset)
                output.write(chunk)
                offset += len(chunk)
                bytes_written += len(chunk)
                if progress_increment is not None:
                    progress_increment(len(chunk))
    return bytes_written


def _progress_accumulator(
    progress: ProgressCallback | None,
) -> Callable[[int], None] | None:
    if progress is None:
        return None

    lock = Lock()
    state = {"bytes_written": 0}

    def increment(chunk_size: int) -> None:
        with lock:
            state["bytes_written"] += chunk_size
            progress(state["bytes_written"], None)

    return increment


def _byte_ranges(content_length: int, threads: int) -> list[tuple[int, int]]:
    segment_size = (content_length + threads - 1) // threads
    ranges = []
    start = 0
    while start < content_length:
        end = min(start + segment_size - 1, content_length - 1)
        ranges.append((start, end))
        start = end + 1
    return ranges


def _throttled_progress(
    progress: ProgressCallback | None,
    content_length: int | None,
    interval_seconds: float,
    clock: Clock,
    last_progress_at: float,
) -> ProgressCallback | None:
    if progress is None:
        return None

    state = {"last_progress_at": last_progress_at}

    def report(bytes_written: int, total: int | None) -> None:
        now = clock()
        effective_total = content_length if content_length is not None else total
        is_complete = effective_total is not None and bytes_written >= effective_total
        if now - state["last_progress_at"] < interval_seconds and not is_complete:
            return
        progress(bytes_written, effective_total)
        state["last_progress_at"] = now

    return report
