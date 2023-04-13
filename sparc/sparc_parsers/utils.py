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


def bisect_and_strip(text, delimiter):
    """split string in 2 at first occurence of a character and remove whitespace
    useful for separating comments from data, keys from values, etc.
    """
    # wrap around to len(text) if not found (-1)
    index = text.find(delimiter) % (len(text) + 1)
    return text[:index].strip(), text[index + len(delimiter) :].strip()
