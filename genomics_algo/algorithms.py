def get_occurences_with_naive_match(pattern, text):
    occurences = []
    if len(pattern) > len(text):
        return occurences
    else:
        for index in range(len(text) - len(pattern) + 1):
            match = True
            for offset in range(len(pattern)):
                if pattern[offset] != text[index + offset]:
                    match = False
                    break
            if match:
                occurences.append(index)
        return occurences
