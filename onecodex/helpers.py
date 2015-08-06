from __future__ import print_function
import sys


SUPPORTED_EXTENSIONS = ["fa", "fasta", "fq", "fastq",
                        "fa.gz", "fasta.gz", "fq.gz", "fastq.gz",
                        "fa.gzip", "fasta.gzip", "fq.gzip", "fastq.gzip"]


def check_for_allowed_file(f):
    for ext in SUPPORTED_EXTENSIONS:
        if f.endswith(ext):
            return True
    stderr("Failed upload: Not an allowed file extension.")
    sys.exit(1)


def stderr(*objs):
    print(*objs, file=sys.stderr)


def is_insecure_platform():
    v = sys.version_info
    if v.major == 2 and v.minor == 7 and v.micro >= 9:
        return False  # >= 2.7.9 includes the new SSL updates
    try:
        import OpenSSL  # noqa
        import ndg  # noqa
        import pyasn1  # noqa
    except ImportError:
        pass
    return True


def warn_if_insecure_platform():
    m = ("\n"
         "######################################################################################\n"
         "#                                                                                    #\n"
         "#  Your version of Python appears to be out of date and lack important security      #\n"
         "#  features. Please update to Python >= 2.7.9 or `pip install requests[security]`.   #\n"
         "#                                                                                    #\n"
         "#  InsecurePlatformWarning: A true SSLContext object is not available. This          #\n"
         "#  prevents urllib3 from configuring SSL appropriately and may cause certain         #\n"
         "#  SSL connections to fail. For more information, see                                #\n"
         "#  https://urllib3.readthedocs.org/en/latest/security.html#insecureplatformwarning.  #\n"
         "#                                                                                    #\n"
         "######################################################################################\n")
    if is_insecure_platform():
        stderr(m)
