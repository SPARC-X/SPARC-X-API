def get_label(fileobj, ext):
    """Return the label of file by stripping the extension (e.g. .ion)"""
    return fileobj.name.rsplit(ext, 1)[0]


def strip_comments(rawtext):
    """ """
    stripped = []
    comments = []
    for line in rawtext.splitlines():
        data, comment = bisect_and_strip(line, "#")
        if data:
            stripped.append(data)
        if comment:
            comments.append(comment)
    return stripped, comments
