import re

__author__ = 'Rusty Gentile'


def code_gen(n_passes, n_elements_per_pass, element_name='e',
             mixer_name='mix', separator_name='sep', model_name='m',
             flowsheet_name='fs', shell_prop_name='prop_fluegas',
             pyomo_environ_name='pe', mixing=True):
    """

    "Programs that write programs are the happiest programs in the world."
    ~~Andrew Hume

    Returns
    -------
    setup_code_block, arc_code_block
    """
    
    total_elements = n_passes * n_elements_per_pass
    setup_code_block = ''
    m, fs, e, mix, sep, pe = model_name, flowsheet_name, element_name, \
                             mixer_name, separator_name, pyomo_environ_name

    for i in range(total_elements):
        setup_code_block += f'{m}.{fs}.{e}{i} = '
        setup_code_block += f'{m}.{fs}.{e}s[{i}]\n'

    setup_code_block += '\n'
    all_elements_line = 'all_elements = ['
    for i in range(total_elements - 1):
        all_elements_line += f'{m}.{fs}.{e}{i}, '

    all_elements_line += f'{m}.{fs}.{e}{total_elements - 1}]'

    def auto_format(line):
        """Breaks up a long line to be PEP8-friendly"""
        res = ''
        while len(line) > 79:
            commas = [l.start() for l in re.finditer(',', line)]
            max_lt_78 = 0
            for c in commas:
                max_lt_78 = c if max_lt_78 < c < 78 else max_lt_78

            if max_lt_78 == 0:
                raise ValueError('Unable to generate formatted code. '
                                 'Most likely, the names provided were too long.')

            res += line[:max_lt_78 + 1] + '\n'
            line = ' ' * 4 + line[max_lt_78 + 2:]

        res += line
        return res

    setup_code_block += auto_format(all_elements_line) + '\n\n'

    shell_in_idx = list(range(total_elements))[-n_elements_per_pass:]
    shell_in_line = 'shell_in_elements = ['
    for i in shell_in_idx[:-1]:
        shell_in_line += f'{m}.{fs}.{e}{i}, '

    shell_in_line += f'{m}.{fs}.{e}{total_elements - 1}]'
    setup_code_block += auto_format(shell_in_line) + '\n\n'
    setup_code_block += f'tube_in_element = {m}.{fs}.{e}0\n'
    setup_code_block += f'tube_out_element = {m}.{fs}.{e}{total_elements-1}\n'

    arc_code_block = '# Add tube-side Arcs\n'
    for i in range(total_elements - 1):
        arc_code_block += f'{m}.{fs}.{e}_tube_Arc{i} = '\
                          f'Arc(source={m}.{fs}.{e}{i}.outlet_1, '\
                          f'destination={m}.{fs}.{e}{i+1}.inlet_1)\n'

    if mixing:
        arc_code_block += '\n# Shell mixer & separator setup\n'
        for i in range(n_passes - 1):
            arc_code_block += f'{m}.{fs}.{mix}{i} = Mixer(default=' \
                              f'{{"property_package": {m}.{fs}.{shell_prop_name},' \
                              f' "num_inlets": {n_elements_per_pass}}})\n'

        arc_code_block += '\n'
        for i in range(n_passes - 1):
            arc_code_block += f'{m}.{fs}.{sep}{i} = Separator(default=' \
                              f'{{"property_package": {m}.{fs}.{shell_prop_name},' \
                              f' "num_outlets": {n_elements_per_pass}}})\n'

        for i in range(n_passes - 1):
            in_start = (n_passes - i - 1) * n_elements_per_pass
            in_end = in_start + n_elements_per_pass - 1
            out_end = in_start - 1
            out_start = in_start - n_elements_per_pass
            inlets = list(range(in_start, in_end + 1))
            outlets = list(range(out_start, out_end + 1))

            arc_code_block += f'\n# Mixer Arcs for pass {i + 1}\n'
            for idx, iidx in enumerate(inlets):
                arc_code_block += f'{m}.{fs}.{e}_shell_pass{i + 1}_m{idx} = Arc(' \
                                  f'source={m}.{fs}.{e}{iidx}.outlet_2, ' \
                                  f'destination={m}.{fs}.{mix}{i}.inlet_{idx + 1})\n'

            arc_code_block += f'\n# Mixer / Separator Arc for pass {i + 1}\n' \
                              f'm.{fs}.{e}_shell_pass{i + 1}_ms_Arc = Arc(' \
                              f'source={m}.{fs}.{mix}{i}.outlet, destination=' \
                              f'{m}.{fs}.{sep}{i}.inlet)\n'

            arc_code_block += f'\n# Separator Arcs for pass {i + 1}\n'
            for idx, iidx in enumerate(outlets):
                arc_code_block += f'{m}.{fs}.{e}_shell_pass{i + 1}_s{idx} = Arc(' \
                                  f'source={m}.{fs}.{sep}{i}.outlet_{idx + 1}, ' \
                                  f'destination={m}.{fs}.{e}{iidx}.inlet_2)\n'

        arc_code_block += f'\n# Apply Arc constraints\n{pe}.TransformationFactory' \
                          f'("network.expand_arcs").apply_to({m}.{fs})\n'

        arc_code_block += '\n# Initialize Mixers and Separators\n'

        for i in range(n_passes - 1):
            arc_code_block += f'{m}.{fs}.{mix}{i}.initialize()\n'

        arc_code_block += '\n'
        for i in range(n_passes - 1):
            arc_code_block += f'{m}.{fs}.{sep}{i}.initialize()\n'

    else:
        for i in range(n_passes - 1):
            in_start = (n_passes - i - 1) * n_elements_per_pass
            in_end = in_start + n_elements_per_pass - 1
            out_end = in_start - 1
            out_start = in_start - n_elements_per_pass
            inlets = list(range(in_start, in_end + 1))
            outlets = list(range(out_start, out_end + 1))
            outlets.reverse()

            arc_code_block += f'\n# Shell Arcs for pass {i + 1}\n'
            for idx in range(len(inlets)):
                ii = inlets[idx]
                oi = outlets[idx]
                arc_code_block += f'{m}.{fs}.{e}_shell_pass{i + 1}_m{idx} = Arc(' \
                                  f'source={m}.{fs}.{e}{ii}.outlet_2, ' \
                                  f'destination={m}.{fs}.{e}{oi}.inlet_2)\n'

        arc_code_block += f'\n# Apply Arc constraints\n{pe}.TransformationFactory' \
                          f'("network.expand_arcs").apply_to({m}.{fs})\n'

    return setup_code_block, arc_code_block


if __name__ == '__main__':
    code = code_gen(8, 7)
    print(code[0])
    print('')
    print(code[1])
