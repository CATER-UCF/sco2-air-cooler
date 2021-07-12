import pyomo.environ as pe


def print_results_0d(e, t=0):
    t_tube_in = pe.value(e.tube.properties_in[t].temperature)
    p_tube_in = pe.value(e.tube.properties_in[t].pressure)
    h_tube_in = pe.value(e.tube.properties_in[t].enth_mol) / 18.01528

    t_tube_out = pe.value(e.tube.properties_out[t].temperature)
    p_tube_out = pe.value(e.tube.properties_out[t].pressure)
    h_tube_out = pe.value(e.tube.properties_out[t].enth_mol) / 18.01528

    t_shell_in = pe.value(e.shell.properties_in[t].temperature)
    p_shell_in = pe.value(e.shell.properties_in[t].pressure)

    t_shell_out = pe.value(e.shell.properties_out[t].temperature)
    p_shell_out = pe.value(e.shell.properties_out[t].pressure)

    heat_duty = pe.value(e.heat_duty[t])
    UA = pe.value(e.overall_heat_transfer_coefficient[t])

    print('--------Tube-----------'
          '')
    print('')
    print(f'T in: {t_tube_in}')
    print(f'T out: {t_tube_out}')
    print('')
    print(f'P in: {p_tube_in}')
    print(f'P out: {p_tube_out}')
    print('')
    print(f'H in: {h_tube_in}')
    print(f'H out: {h_tube_out}')
    print('')
    print('--------Shell-----------')
    print('')
    print(f'T in: {t_shell_in}')
    print(f'T out: {t_shell_out}')
    print('')
    print(f'P in: {p_shell_in}')
    print(f'P out: {p_shell_out}')
    print('')
    print(f'Overall HTC: {UA}')
    print(f'Heat duty (MW): {heat_duty * 1e-6}')

    return heat_duty * 1e-6
