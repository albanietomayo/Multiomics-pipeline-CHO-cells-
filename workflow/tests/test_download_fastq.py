#!/usr/bin/env python3
import hashlib
import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflow/scripts/download_fastq.py"


def load():
    sys.modules.setdefault("pandas", types.SimpleNamespace())
    spec = importlib.util.spec_from_file_location("download_fastq_test_support", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DOWNLOAD = load()


class Response:
    def __init__(self, chunks, status_code=200, headers=None):
        self.chunks = chunks
        self.status_code = status_code
        self.headers = headers or {}

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        del chunk_size
        for chunk in self.chunks:
            if isinstance(chunk, Exception):
                raise chunk
            yield chunk


class DownloadFastqTests(unittest.TestCase):
    @mock.patch.object(DOWNLOAD.time, "sleep")
    def test_interrupted_https_download_resumes_and_validates(self, _sleep):
        payload = b"abcdefghij"
        responses = [
            Response([payload[:4], OSError("truncated")]),
            Response(
                [payload[4:]],
                status_code=206,
                headers={"Content-Range": "bytes 4-9/10"},
            ),
        ]

        with tempfile.TemporaryDirectory() as temporary, mock.patch.object(
            DOWNLOAD.requests, "get", side_effect=responses
        ) as get:
            output = Path(temporary) / "reads.fastq.gz"
            DOWNLOAD.download_and_validate(
                "ftp.sra.ebi.ac.uk/reads.fastq.gz",
                output,
                len(payload),
                hashlib.md5(payload).hexdigest(),
                retry_delay=0,
            )

            self.assertEqual(output.read_bytes(), payload)
            self.assertFalse(Path(str(output) + ".part").exists())
            self.assertEqual(get.call_args_list[1].kwargs["headers"], {"Range": "bytes=4-"})

    @mock.patch.object(DOWNLOAD.time, "sleep")
    def test_server_ignoring_range_restarts_without_duplication(self, _sleep):
        payload = b"abcdefghij"
        responses = [
            Response([payload[:4], OSError("truncated")]),
            Response([payload]),
        ]

        with tempfile.TemporaryDirectory() as temporary, mock.patch.object(
            DOWNLOAD.requests, "get", side_effect=responses
        ):
            output = Path(temporary) / "reads.fastq.gz"
            DOWNLOAD.download_and_validate(
                "https://ftp.sra.ebi.ac.uk/reads.fastq.gz",
                output,
                len(payload),
                hashlib.md5(payload).hexdigest(),
                retry_delay=0,
            )

            self.assertEqual(output.read_bytes(), payload)


if __name__ == "__main__":
    unittest.main()
