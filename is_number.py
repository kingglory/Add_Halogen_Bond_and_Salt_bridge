# write a function to make sure the input is number (to choose the right line that we will use to compute the distance and angles

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False

# test is_number()

print is_number(2)
