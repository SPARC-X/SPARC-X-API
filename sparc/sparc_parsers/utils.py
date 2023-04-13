from warnings import warn

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

def read_block_input(block, validator=None):
    """Read blocks of inputs from ion or inpt file and convert with validator
    """
    block_dict = {}
    multiline_key = ""
    use_validator = True if validator else False
    for line in block:
        if ":" not in line:
            # no key, assume multiline value
            block_dict[multiline_key].append(line.strip())
            continue
        key, value = bisect_and_strip(line, ":")
        key = key.upper()
        # print(key, value)
        if key and value:
            block_dict[key] = value
        elif key:
            # no value, assume that this key has a list of values
            # in the following lines
            block_dict[key] = []
            multiline_key = key
    for key, val in block_dict.items():
        # print(key, val)
        _use_validator_this_key = use_validator
        if _use_validator_this_key:
            if key not in validator.parameters.keys():
                warn(f"Key {key} not in validator's parameter list, ignore value conversion!")
                _use_validator_this_key = False
        if _use_validator_this_key:
            val = validator.convert_string_to_value(key, val)
        block_dict[key] = val
    return block_dict
