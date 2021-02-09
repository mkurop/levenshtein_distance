# -*-coding; utf-8 -*-

"""Module for computation of the levenshtein distance between two strings"""

__author__ = 'markurop'

import functools


@functools.total_ordering
class hpt:

    def __eq__(self, other):
        return ((self.pos_in, self.pos_true, self.cost) == (other.pos_in, other.pos_true, other.cost))

    def __le__(self, other):
        return ((self.pos_in, self.pos_true, self.cost) < (other.pos_in, other.pos_true, other.cost))

    def correct(self,h,bck_ptr):
        self.pos_true = h.pos_true+1
        self.pos_in = h.pos_in+1
        self.cost = h.cost
        self.back_ptr = bck_ptr
        self.num_del = h.num_del
        self.num_ins = h.num_ins
        self.num_sub = h.num_sub

    def substitution(self,h,bck_ptr):
        self.pos_true = h.pos_true+1
        self.pos_in = h.pos_in+1
        self.cost = h.cost+1
        self.back_ptr = bck_ptr
        self.num_del = h.num_del
        self.num_ins = h.num_ins
        self.num_sub = h.num_sub+1

    def insertion(self,h,bck_ptr):
        self.pos_true = h.pos_true
        self.pos_in = h.pos_in+1
        self.cost = h.cost+1
        self.back_ptr = bck_ptr
        self.num_del = h.num_del
        self.num_ins = h.num_ins+1
        self.num_sub = h.num_sub

    def deletion(self,h,bck_ptr):
        self.pos_true = h.pos_true+1
        self.pos_in = h.pos_in
        self.cost = h.cost+1
        self.back_ptr = bck_ptr
        self.num_del = h.num_del+1
        self.num_ins = h.num_ins
        self.num_sub = h.num_sub

    def start_correct(self):
        self.pos_true = 1
        self.pos_in = 1
        self.cost = 0
        self.back_ptr = None
        self.num_del = 0
        self.num_ins = 0
        self.num_sub = 0

    def start_substitution(self):
        self.pos_true = 1
        self.pos_in = 1
        self.cost = 1
        self.back_ptr = None
        self.num_del = 0
        self.num_ins = 0
        self.num_sub = 1

    def start_insertion(self):
        self.pos_true = 0
        self.pos_in = 1
        self.cost = 1
        self.back_ptr = None
        self.num_del = 0
        self.num_ins = 1
        self.num_sub = 0

    def start_deletion(self):
        self.pos_true = 1
        self.pos_in = 0
        self.cost = 1
        self.back_ptr = None
        self.num_del = 1
        self.num_ins = 0
        self.num_sub = 0

    def __init__(self,error_type,bck_ptr,prev_hpt,start_flag):
        """
        :param error_type: 'correct','substitution','deletion','insertion'
        :param bck_ptr: back pointer, None for (start_flag == True)
        :param prev_hpt: previous hpt, None for (start_flag == True)
        :param start_flag: True for first invocation False otherwise
        :return:
        """
        if start_flag:
            start_dict = {'correct' : self.start_correct,
                         'substitution' : self.start_substitution,
                         'insertion' : self.start_insertion,
                         'deletion' : self.start_deletion}
            start_dict[error_type]()
        else:
            non_start_dict = {'correct' : self.correct,
                         'substitution' : self.substitution,
                         'insertion' : self.insertion,
                         'deletion' : self.deletion}
            non_start_dict[error_type](prev_hpt,bck_ptr)

def init_hpt(true_itm,in_itm):
    hpt_lst = []
    if true_itm != in_itm:
        hpt_lst.append(hpt('substitution',None,None,True))
        hpt_lst.append(hpt('insertion',None,None,True))
        hpt_lst.append(hpt('deletion',None,None,True))
    else:
        hpt_lst.append(hpt('correct',None,None,True))

    return hpt_lst

def expand(hpt_lst,true_seq,in_seq):
    hpt_lst_out = []
    no_more_to_expand = True
    for bck_ptr,h in enumerate(hpt_lst):
        if h.pos_true >= len(true_seq) and h.pos_in < len(in_seq):
            hpt_lst_out.append(hpt('insertion',bck_ptr,h,False))
            no_more_to_expand = False
        if h.pos_true < len(true_seq) and h.pos_in >= len(in_seq):
            hpt_lst_out.append(hpt('deletion',bck_ptr,h,False))
            no_more_to_expand = False
        if h.pos_true < len(true_seq) and h.pos_in < len(in_seq):
            if true_seq[h.pos_true] != in_seq[h.pos_in]:
                hpt_lst_out.append(hpt('substitution',bck_ptr,h,False))
                hpt_lst_out.append(hpt('insertion',bck_ptr,h,False))
                hpt_lst_out.append(hpt('deletion',bck_ptr,h,False))
            else:
                hpt_lst_out.append(hpt('correct',bck_ptr,h,False))
            no_more_to_expand = False
        if h.pos_true >= len(true_seq) and h.pos_in >= len(in_seq):
            hpt_lst_out.append(h)
    return no_more_to_expand, hpt_lst_out

def recombination(hpt_lst):
    #single best
    hpt_lst_sorted = sorted(hpt_lst)
    hpt_lst_out = [hpt_lst_sorted[0]]
    for i in range(len(hpt_lst_sorted)-1):
        if ((hpt_lst_sorted[i].pos_in, hpt_lst_sorted[i].pos_true) != (hpt_lst_sorted[i+1].pos_in, hpt_lst_sorted[i+1].pos_true)):
            hpt_lst_out.append(hpt_lst_sorted[i+1])
    return hpt_lst_out

class levenshtein:

    """Class for computing Levenshtein distance between two strings"""

    def __init__(self):

        pass

    def compare(self,true_seq,input_seq):

        """The actual comparing method

        :param true_seq: true string
        :param input_seq: the string to be compared with true string

        :returns:
            :num_sub: number of substitutions
            :num_del: number of deletions
            :num_ins: number of insertions

        """

        # start program here
        hpt_lst = init_hpt(true_seq[0], input_seq[0])

        len_true = len(true_seq)
        len_in = len(input_seq)
        no_more_to_expand = False
        seq_hpt_lst = []
        while not no_more_to_expand:
            no_more_to_expand, hpt_lst = expand(hpt_lst, true_seq, input_seq)
            hpt_lst = recombination(hpt_lst)
            if not no_more_to_expand:
                seq_hpt_lst.append(hpt_lst)

        # backtracking


        return seq_hpt_lst[-1][0].num_sub, seq_hpt_lst[-1][0].num_del, seq_hpt_lst[-1][0].num_ins

if __name__ == "__main__":

    string1 = ['I','am','a','glorious','cat']
    string2 = ['I','am','a','dog']

    lev = levenshtein()

    num_sub, num_del, num_ins = lev.compare(string1,string2)

    print(f"Number of substitutions: {num_sub}, number of deletions: {num_del}, number of insertions: {num_ins}")


