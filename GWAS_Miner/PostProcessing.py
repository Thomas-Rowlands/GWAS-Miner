

def clean_output_annotations(study: object) -> object:
    """
    Trim unwanted characters from study annotations.
    """
    for passage in study['documents'][0]['passages']:
        for annot in passage["annotations"]:
            start_trimmed, end_trimmed = 0, 0
            cleaned = False
            is_p_value = annot.infons["type"] == "significance"
            while not cleaned:
                if not annot.text[-1].isdigit() and not annot.text[-1].isalpha():
                    annot.text = annot.text[:-1]
                    end_trimmed += 1
                elif (not annot.text[0].isdigit() and not annot.text[0].isalpha()) or (not annot.text[0].isdigit() and is_p_value):
                    annot.text = annot.text[1:]
                    start_trimmed += 1
                else:
                    cleaned = True
                    annot.locations[0].offset += start_trimmed
                    annot.locations[0].length = annot.locations[0].length - (start_trimmed + end_trimmed)
    return study
