import re
import numpy as np

def tokenizer(smile):
    "Tokenizes SMILES string"
    pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|_|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regezz = re.compile(pattern)
    tokens = [token for token in regezz.findall(smile)]
    assert smile == ''.join(tokens), ("{} could not be joined".format(smile))
    return tokens

def chiralizer(smile):
    toks = tokenizer(smile)
    beta_atom = toks[1]
    chiral_center = True
    if beta_atom == 'C' or beta_atom == 'Si':
        pass
    else:
        chiral_center = False

    if chiral_center:
        try:
            branch1_start = toks[2]
            if branch1_start == '(':
                branch_idx = 3
                additional_branch = False
                branch1_ended = False
                branch1 = []

                while not branch1_ended:
                    next_tok = toks[branch_idx]
                    if next_tok == '(':
                        additional_branch = True
                    elif next_tok == ')' and additional_branch:
                        additional_branch = False
                    elif next_tok == ')' and not additional_branch:
                        branch1_ended = True
                    else:
                        branch1.append(next_tok)
                    branch_idx += 1

                branch2 = toks[branch_idx:]

                if branch1 == branch2:
                    chiral_center = False
                else:
                    toks[1] = '[{}@H]'.format(beta_atom)
                    enantiomer1 = ''.join(toks)
                    toks[1] = '[{}@@H]'.format(beta_atom)
                    enantiomer2 = ''.join(toks)

            else:
                chiral_center = False

        except IndexError:
            chiral_center = False



    if chiral_center:
        return chiral_center, [enantiomer1, enantiomer2]
    else:
        return chiral_center, ''.join(toks)
