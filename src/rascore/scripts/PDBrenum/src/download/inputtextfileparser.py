"""input text file will be parsed by comma, space, tab and new line"""


def input_text_file_parser(filename):
    with open(filename) as file:
        file_contents = [line.rstrip() for line in file]
    comma_space_tab_newline_parser = list()
    for n_element in file_contents:
        if "," in n_element:
            n_parsed = n_element.split(",")
            comma_space_tab_newline_parser.append(n_parsed)
        elif " " in n_element:
            n_parsed = n_element.split(" ")
            comma_space_tab_newline_parser.append(n_parsed)
        elif "\t" in n_element:
            n_parsed = n_element.split("\t")
            comma_space_tab_newline_parser.append(n_parsed)
        elif "\n" in n_element:
            n_parsed = n_element.split("\n")
            comma_space_tab_newline_parser.append(n_parsed)
        else:
            comma_space_tab_newline_parser.append(n_element)

    if_list_parse_list = list()
    for n_element in comma_space_tab_newline_parser:
        if type(n_element) == list:
            for d_list_parsed in n_element:
                if d_list_parsed not in if_list_parse_list:
                    if_list_parse_list.append(d_list_parsed)

    final_parsed_list = list()
    for n_elem in if_list_parse_list:
        for _ in range(3):
            n_elem = n_elem.strip(",")
            n_elem = n_elem.strip()

        if len(n_elem) < 4:
            continue
        final_parsed_list.append(n_elem)

    return final_parsed_list
