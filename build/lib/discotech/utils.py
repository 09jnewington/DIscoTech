import logging as lg
import sys


def setup_logger(
    logfile: str, file_level: int = lg.DEBUG, stream_level: int = lg.INFO
) -> lg.Logger:
    """
    Sets up a generic logger that logs `file_level` messages to the file specified by `logfile` and `stream_level`
    messages to standard output.

    :param logfile: Filename for the log file.
    :param file_level: Level of messages logged to the log file.
    :param stream_level: Level of messages logged to stdout.
    :return: Logger.
    """

    logger = lg.getLogger(__name__)

    formatter = lg.Formatter("%(asctime)s %(levelname)-8s: %(message)s")
    filehandler = lg.FileHandler(logfile)
    filehandler.setLevel(file_level)
    filehandler.setFormatter(formatter)
    streamhandler = lg.StreamHandler(sys.stdout)
    streamhandler.setLevel(stream_level)
    streamhandler.setFormatter(formatter)
    logger.addHandler(filehandler)
    logger.addHandler(streamhandler)

    logger_level = file_level if file_level <= stream_level else stream_level

    logger.setLevel(logger_level)
    logger.propagate = False

    return logger
