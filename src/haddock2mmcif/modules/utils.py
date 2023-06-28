def backmap(dict, value):
    """Find a key in a dictionary by its value."""
    return list(dict.keys())[list(dict.values()).index(value)]
