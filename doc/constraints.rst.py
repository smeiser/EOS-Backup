import eos
import re
import yaml
from jinja_util import print_template

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

qn_to_link_map = {
    ord(':'): 'co', ord('@'): 'at', ord('/'): 'sl', ord('_'): 'un',
    ord('('): 'po', ord(')'): 'pc', ord('+'): 'pp', ord('-'): 'mm',
    ord('>'): 'to'
}

def make_constraints():
    result = []
    for qn, entry in eos.Constraints():
        data = yaml.safe_load(entry.serialize())

        if 'observable' in data:
            unique_observables = { str(data['observable']) }
        elif 'observables' in data:
            unique_observables = { str(o) for o in data['observables'] }

        observable_entries = []
        for oqn in unique_observables:
            anchor = oqn.translate(qn_to_link_map).lower()
            observable_entries.append(f'`{ oqn } <observables.html#{ anchor }>`_')

        constraint = {
            'observables': ', '.join(observable_entries)
        }
        result.append((qn, constraint))

    return result

if __name__ == '__main__':

    print_template(__file__,
        version = eos.__version__,
        constraints = make_constraints(),
        len = len,
    )
