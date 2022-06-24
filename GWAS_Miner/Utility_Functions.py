from spacy.tokens import Span


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
    def remove_duplicate_associations(assoc_input):
        filtered_results = []
        for i in assoc_input:
            is_unique = True
            marker_idx = i.marker.token.idx if type(i.marker.token) != Span else i.marker.token.ents[0].start
            significance_idx = i.significance.token.idx if type(i.significance.token) != Span else i.significance.token.ents[0].start
            phenotype_idx = i.phenotype.token.idx if type(i.phenotype.token) != Span else i.phenotype.token.ents[0].start
            for l in filtered_results:
                if i == l:
                    is_unique = False
                    break
                used_marker_idx = l.marker.token.idx if type(l.marker.token) != Span else l.marker.token.ents[0].start
                used_significance_idx = l.significance.token.idx if type(l.significance.token) != Span else l.significance.token.ents[0].start
                used_phenotype_idx = l.phenotype.token.idx if type(l.phenotype.token) != Span else l.phenotype.token.ents[0].start
                if marker_idx == used_marker_idx and significance_idx == used_significance_idx and phenotype_idx == used_phenotype_idx:
                    is_unique = False
                    break
            if is_unique:
                filtered_results.append(i)
        return filtered_results

    @staticmethod
    def remove_duplicate_bioc_associations(assoc_input: list) -> list:
        filtered_results = []
        for i in assoc_input:
            is_unique = True
            for l in filtered_results:
                if i == l:
                    is_unique = False
                    break
                if i.nodes == l.nodes:
                    is_unique = False
                    break
            if is_unique:
                filtered_results.append(i)
        return filtered_results

    @staticmethod
    def retrieve_value_indexes(search, list_input):
        output = []
        for item in list_input:
            if item == search:
                output.append(list_input.index(search))
                list_input.remove(search)
        return output

    @staticmethod
    def expand_xpath_output(xpath_result):
        if type(xpath_result) == list:
            if xpath_result == []:
                return None
            else:
                return xpath_result[0]
        else:
            return xpath_result
