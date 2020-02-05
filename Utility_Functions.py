

class Utility:
    @staticmethod
    def string_array_compare(array_a, array_b):
        for entry in array_a:
            if entry not in array_b:
                return False
        return True

    @staticmethod
    def remove_duplicates(list_input):
        output = []
        for item in list_input:
            if item not in output:
                output.append(item)
        return output

    @staticmethod
    def retrieve_value_indexes(search, list_input):
        output = []
        for item in list_input:
            if item == search:
                output.append(list_input.index(search))
                list_input.remove(search)
        return output

