"""
Mixins for extended functionality in the new model system.
"""

import os
import os.path
from typing import Optional, BinaryIO, Union

import click
import httpx

from onecodex.exceptions import OneCodexException


class DownloadMixin:
    """Mixin for models that support downloading files."""
    
    def download(
        self, 
        path: Optional[str] = None, 
        file_obj: Optional[BinaryIO] = None, 
        progressbar: bool = False
    ) -> Optional[str]:
        """Download files from One Codex.
        
        Args:
            path: Full path to save the file to. If omitted, defaults to the original filename
                 in the current working directory.
            file_obj: Rather than save the file to a path, write it to this file-like object.
            progressbar: Display a progress bar using Click for the download?
            
        Returns:
            The path the file was downloaded to, if applicable. Otherwise, None.
        """
        return self._download(
            "download_uri",
            getattr(self, 'filename', None),
            path=path,
            file_obj=file_obj,
            progressbar=progressbar,
        )
    
    def _download(
        self,
        endpoint_method: str,
        filename: Optional[str] = None,
        path: Optional[str] = None,
        file_obj: Optional[BinaryIO] = None,
        progressbar: bool = False,
    ) -> Optional[str]:
        """Internal download implementation."""
        
        # Check model state
        if hasattr(self, 'visibility') and getattr(self, 'visibility') == "awaiting data":
            raise OneCodexException("Sample has not finished processing. Please try again later.")
        
        if path and file_obj:
            raise OneCodexException("Please specify only one of: path, file_obj")
        
        # Get download URL from API
        self._check_bound()
        endpoint = f"{self._resource_path}/{self.id}/{endpoint_method}"
        
        try:
            download_info = self._client.http_client.post(endpoint)
        except Exception as e:
            raise OneCodexException(f"Failed to get download URL: {e}")
        
        # Determine output path
        if path is None and file_obj is None:
            if filename is None:
                if "save_as_filename" not in download_info:
                    raise OneCodexException("Please specify `path`, `file_obj`, or provide a filename.")
                filename = download_info["save_as_filename"]
            path = os.path.join(os.getcwd(), filename)
        
        if path and os.path.exists(path):
            raise OneCodexException(f"{path} already exists. Will not overwrite.")
        
        # Download the file
        download_url = download_info["download_uri"]
        
        try:
            # Create a new httpx client for the download with retry logic
            with httpx.Client(
                timeout=300.0,  # 5 minute timeout for downloads
                transport=httpx.HTTPTransport(retries=5)
            ) as client:
                
                with client.stream("GET", download_url) as response:
                    response.raise_for_status()
                    
                    # Open output file
                    if path:
                        f_out = open(path, "wb")
                    else:
                        f_out = file_obj
                    
                    try:
                        if progressbar:
                            # Get file size from response headers
                            total_size = int(response.headers.get('content-length', 0))
                            progress_label = os.path.basename(path) if path else filename or "download"
                            
                            with click.progressbar(length=total_size, label=progress_label) as bar:
                                for chunk in response.iter_bytes(chunk_size=8192):
                                    bar.update(len(chunk))
                                    f_out.write(chunk)
                        else:
                            for chunk in response.iter_bytes(chunk_size=8192):
                                f_out.write(chunk)
                    finally:
                        # Only close if we opened the file
                        if path:
                            f_out.close()
                        
        except KeyboardInterrupt:
            if path and os.path.exists(path):
                os.remove(path)
            raise
        except httpx.HTTPStatusError as exc:
            if exc.response.status_code == 401:
                raise OneCodexException("You must be logged in to download files.")
            elif exc.response.status_code == 402:
                raise OneCodexException(
                    "You must either have a premium platform account or be in "
                    "a notebook environment to download files. Please contact us "
                    "about your subscription at support@onecodex.com."
                )
            elif exc.response.status_code == 403:
                raise OneCodexException("You are not authorized to download this file.")
            else:
                raise OneCodexException(f"Download failed with HTTP status code {exc.response.status_code}.")
        
        return path


class UploadMixin:
    """Mixin for models that support uploading files."""
    
    @classmethod
    def upload(
        cls, 
        file_path: str, 
        progressbar: Optional[click.progressbar] = None,
        **kwargs
    ):
        """Upload a file to One Codex.
        
        Args:
            file_path: Path to the file to upload
            progressbar: Optional Click progressbar for progress display
            **kwargs: Additional upload parameters
            
        Returns:
            The created model instance
        """
        cls._check_bound()
        
        if not os.path.exists(file_path):
            raise OneCodexException(f"File not found: {file_path}")
        
        # This is a simplified implementation
        # The full implementation would need to handle:
        # - Multipart uploads for large files
        # - Pre-upload validation
        # - Progress tracking
        # - Different upload types (samples vs documents)
        
        # For now, this is a placeholder that shows the structure
        raise NotImplementedError(
            "Upload functionality not yet fully implemented in new model system. "
            "This requires porting the complex multipart upload logic from lib/upload.py"
        )


def truncate_string(s: str, length: int = 24) -> str:
    """Truncate a string to the specified length with ellipsis."""
    if len(s) < length - 3:
        return s
    else:
        s = s[0:(length - 3)]
        if s[-1] == ".":
            s = s + ".."
        else:
            s = s + "..."
        return s