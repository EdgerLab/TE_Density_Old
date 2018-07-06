# By Scott teresi
#------------------------------------------------------------------
from genic_elements import *
from collections import deque, defaultdict
import os
import csv
#------------------------------------------------------------------

def get_densities(genes,transposons,window,increment,max_window,gtf_inputfile):

    #file_list = ['gtf_chunks_aa.gtf','gtf_chunks_ab.gtf','gtf_chunks_ac.gtf','gtf_chunks_ad.gtf','gtf_chunks_ae.gtf','gtf_chunks_af.gtf']
    file_list = ['xx00','xx01','xx02','xx03','xx04','xx05','xx06']

    for item in file_list:
        if gtf_inputfile == item:
            #item = item[:-3]
            overlaps  = 'overlaps_' + item + '.csv'

    try:
        os.remove(overlaps)
    except FileNotFoundError:
        pass

    sp = ','
    f = open(overlaps,'w+')
    print('g_start' + sp + 'g_stop' + sp + 'g_chromosome' + sp + 't_start' + sp + 't_stop' + sp + 't_chromosome' + sp + 't_type' + sp + 't_family',file=f)
    f.close()

    dict_type_list = []
    dict_fam_list = []
    for elem in transposons:
        if elem.getFamily() not in dict_fam_list:
            dict_fam_list.append(elem.getFamily())
        if elem.getTe_Type() not in dict_type_list:
            dict_type_list.append(elem.getTe_Type())

    while window <= max_window:
        for elem in genes:
            g_start = elem.getStart()
            g_stop = elem.getStop()
            g_length = elem.getLength()
            g_chromosome = elem.getChromosome()
            g_number = elem.getNumber()
            prox_left = elem.getProxLeft()
            prox_right = elem.getProxRight()

            left_window_start = g_start - window
            left_window_stop = g_start - 1
            if left_window_start < 0 :
                left_window_start = 0
            right_window_start= g_stop + 1
            right_window_stop = g_stop + window
            g_adj_length = right_window_stop - left_window_start + 1


            for poson in transposons:
                t_start = poson.getStart()
                t_stop = poson.getStop()
                t_length = poson.getLength()
                t_chromosome = poson.getChromosome()
                t_type = poson.getTe_Type()
                t_family = poson.getFamily()
                data_safety = 0

                left_density = 0
                intra_density = 0
                right_density = 0

                if g_chromosome == t_chromosome and g_stop >= g_start and t_stop >= t_start:
                # the data safety variable is used for simple calculations

                    if t_start <= left_window_start and g_start <= t_stop <= g_stop: # The TE began outside the window, but ended inside the gene
                        left_density = left_outside_from_window(window,left_window_stop,left_window_start,left_density)
                        intra_density = left_inside(t_stop,g_start,g_stop,t_type,intra_density)

                        #elem.prox_left = L_prox_check(left_window_stop,t_stop,prox_left)
                        elem.prox_left = 0 # Made the proximity 0 due to overlap

                        data_safety = 0 # special case, no more writing

                    elif t_start <= left_window_start and left_window_start <= t_stop <= g_start: # TE began outside window but ended inside window, not in gene
                        t_length = t_stop - left_window_start + 1
                        if t_stop == left_window_start: # I am in a situation where the end of the TE is right on the spot where the window is placed. Would return a 0 for length, but will alter to return a 1.
                            t_length = 1
                        elem.prox_left = L_prox_check(left_window_stop,t_stop,prox_left)
                        data_safety = 1 # calculates density inside the window

                    elif left_window_start <= t_start <= g_start and left_window_start <= t_stop <= g_start: # TE is before the gene start and no blending, but within window, purely within the window
                        elem.prox_left = L_prox_check(left_window_stop,t_stop,prox_left)
                        data_safety = 1

                    elif t_start >= left_window_start and t_start <= g_start and t_stop <= g_stop and t_stop >= g_start: # TE is starts before the gene, and ends within the gene, positive for blending,
                        left_density = left_outside_in_window(window,left_window_stop,t_start,left_density)
                        intra_density = left_inside(t_stop,g_start,g_stop,t_type,intra_density)
                        #elem.prox_left = L_prox_check(left_window_stop,t_stop,prox_left)
                        elem.prox_left = 0 # Made the proximity 0 due to overlap
                        data_safety = 0 # special case, no more writing

                    #------------
                    # center

                    elif t_start >= g_start and t_stop <= g_stop: #TE is within the gene
                        # Only need to compute once and it is 1 value per gene
                        elem.they_are_inside += 1 # Add a count to the number of TE's solely within the gene.
                        data_safety = 2 # Flagging for later tagging.

                    #-------------
                    # right side

                    elif g_start <= t_start <= g_stop and g_stop <= t_stop <= right_window_stop: # TE starts inside of the gene, and ends outside, positive for blending, does not extend past window
                        right_density = right_outside_in_window(window,right_window_start,t_stop,right_density)
                        intra_density = right_inside(t_start,g_start,g_stop,t_type,intra_density)
                        #elem.prox_right = R_prox_check(right_window_start,t_start,prox_right)
                        elem.prox_right = 0 # Made proximity 0 due to overlap
                        data_safety = 0  # special case, no more writing

                    elif g_start <= t_start <= g_stop and t_stop >= right_window_stop: # TE starts inside of the gene, and ends past the window
                        right_density = right_outside_from_window(window,right_window_start,right_window_stop,right_density)
                        intra_density = right_inside(t_start,g_start,g_stop,t_type,intra_density)
                        #elem.prox_right = R_prox_check(right_window_start,t_start,prox_right)
                        elem.prox_right = 0 # Made proximity 0 due to overlap
                        data_safety = 0  # special case, no more writing

                    elif g_stop <= t_start <= right_window_stop and t_stop >= right_window_stop: # TE sits upstream of the gene, and extends past the window
                        t_length = right_window_stop - t_start + 1 # Update te_length to reflect the portion of the bases in the window that are actually TE
                        if t_start == right_window_stop: # I am in a situation where the start of the TE is right on the spot where the window is placed. Would return a 0 for length, but will alter to return a 1.
                            t_length = 1
                        elem.prox_right = R_prox_check(right_window_start,t_start,prox_right)
                        data_safety = 3

                    elif g_stop <= t_start <= right_window_stop and g_stop <= t_stop <= right_window_stop:
                        elem.prox_right = R_prox_check(right_window_start,t_start,prox_right)
                        data_safety = 3


                    #-----------------------------------------------------------------------------------------------
                    # exceptional case where a TE covers a gene from left to right
                    elif left_window_start <= t_start <= g_start and g_stop <= t_stop <= right_window_stop: # TE covers from left to right but within the window, exceptional 1 in notebook
                        #left_density = left_outside_in_window(window,left_window_stop,t_start,left_density)
                        #right_density = right_outside_in_window(window,right_window_start,t_stop,right_density)
                        #intra_density += 1 # because the entire inside is covered
                        data_safety = 0

                        f = open(overlaps,'a+')
                        sp = ','
                        print(str(g_start) + sp + str(g_stop) + sp + str(g_chromosome) + sp + str(t_start) + sp + str(t_stop) + sp + str(t_chromosome) + sp + str(t_type) + sp +str(t_family),file=f)
                        f.close()

                    elif t_start <= left_window_start and g_stop <= t_stop <= right_window_stop: # TE comes from outside the window and overlaps the gene, exceptional #2
                        #left_density= left_outside_from_window(window,left_window_stop,left_window_start,left_density)
                        #right_density= right_outside_in_window(window,right_window_start,t_stop,right_density)
                        #intra_density += 1 # because the entire inside is covered
                        data_safety = 0

                        f = open(overlaps,'a+')
                        sp = ','
                        print(str(g_start) + sp + str(g_stop) + sp + str(g_chromosome) + sp + str(t_start) + sp + str(t_stop) + sp + str(t_chromosome) + sp + str(t_type) + sp +str(t_family),file=f)
                        f.close()

                    elif left_window_start <= t_start <= g_start and right_window_stop <= t_stop:
                        #left_density = left_outside_in_window(window,left_window_stop,t_start,left_density)
                        #right_density = right_outside_from_window(window,right_window_start,right_window_stop,right_density)
                        #intra_density += 1 # because the entire inside is covered
                        data_safety = 0

                        f = open(overlaps,'a+')
                        sp = ','
                        print(str(g_start) + sp + str(g_stop) + sp + str(g_chromosome) + sp + str(t_start) + sp + str(t_stop) + sp + str(t_chromosome) + sp + str(t_type) + sp +str(t_family),file=f)
                        f.close()

                    else:
                        # This condition happens when there is no match and we don't want to add anything
                        data_safety = 0
                    #-----------------------------------------------------------------------------------------------
                    if data_safety == 1:
                        # Left aggregate, here I will add to my aggregate because it passed the check.
                        left_density = t_length/window

                    elif data_safety == 3:
                        # Right aggregate, here I will add to my aggregate because it passed the check.
                        right_density = t_length/window

                    elif data_safety == 2:
                        intra_density = t_length/(g_stop-g_start+1)

#-------------------------------------------------------------------------------------


                    # Here I will update the densities
                    if t_type == 'DNA':
                        elem.DNA_left += left_density
                        elem.DNA_intra += intra_density
                        elem.DNA_right += right_density
                    elif t_type == 'LTR':
                        elem.LTR_left += left_density
                        elem.LTR_intra += intra_density
                        elem.LTR_right += right_density
                    elif t_type == 'LINE':
                        elem.LINE_left += left_density
                        elem.LINE_intra += intra_density
                        elem.LINE_right += right_density
                    elif t_type == 'Unknown':
                        elem.Unknown_left += left_density
                        elem.Unknown_intra += intra_density
                        elem.Unknown_right += right_density

                    # Families

                    if t_family == 'LTRUnknown':
                        elem.LTRUnknown_left += left_density
                        elem.LTRUnknown_intra += intra_density
                        elem.LTRUnknown_right += right_density
                    elif t_family == 'MULE':
                        elem.MULE_left += left_density
                        elem.MULE_intra += intra_density
                        elem.MULE_right += right_density
                    elif t_family == 'Gypsy':
                        elem.Gypsy_left += left_density
                        elem.Gypsy_intra += intra_density
                        elem.Gypsy_right += right_density
                    elif t_family == 'Unknown':
                        elem.Unknown_fam_left += left_density
                        elem.Unknown_fam_intra += intra_density
                        elem.Unknown_fam_right += right_density
                    elif t_family == 'CMC_EnSpm':
                        elem.CMC_EnSpm_left += left_density
                        elem.CMC_EnSpm_intra += intra_density
                        elem.CMC_EnSpm_right += right_density
                    elif t_family == 'Copia':
                        elem.Copia_left += left_density
                        elem.Copia_intra += intra_density
                        elem.Copia_right += right_density
                    elif t_family == 'LINE':
                        elem.LINE_fam_left += left_density
                        elem.LINE_fam_intra += intra_density
                        elem.LINE_fam_right += right_density
                    elif t_family == 'hAT':
                        elem.hAT_left += left_density
                        elem.hAT_intra += intra_density
                        elem.hAT_right += right_density
                    elif t_family == 'PIF_Harbinger':
                        elem.PIF_Harbinger_left += left_density
                        elem.PIF_Harbinger_intra += intra_density
                        elem.PIF_Harbinger_right += right_density

        for item in file_list:
            if gtf_inputfile == item:
                moving_filename = str(window) + '_' + item + '_density_data.csv'

                with open(moving_filename, 'w') as f_out:
                    fieldnames = ['maker_name', 'chromosome', 'start', 'stop', 'length',
                    'DNA_left', 'DNA_intra', 'DNA_right', 'LTR_left', 'LTR_intra', 'LTR_right',
                    'Unknown_left', 'Unknown_intra', 'Unknown_right', 'LINE_left', 'LINE_intra', 'LINE_right',
                    'MULE_left', 'MULE_intra', 'MULE_right', 'Gypsy_left', 'Gypsy_intra', 'Gypsy_right',
                    'Unknown_fam_left','Unknown_fam_intra', 'Unknown_fam_right', 'CMC_EnSpm_left',
                    'CMC_EnSpm_intra', 'CMC_EnSpm_right', 'Copia_left','Copia_intra','Copia_right',
                    'LINE_fam_left','LINE_fam_intra','LINE_fam_right','hAT_left','hAT_intra','hAT_right',
                    'PIF_Harbinger_left','PIF_Harbinger_intra','PIF_Harbinger_right','prox_left','prox_right',
                    'LTRUnknown_left','LTRUnknown_intra','LTRUnknown_right',
                    'they_are_inside','number', 'window_size']# Supply the header to the file, this also orders your output columns in the supplied order. Must have all attributes listed.
                    # You will need to edit the fieldnames to match your data

                    f_out = csv.DictWriter(f_out,fieldnames=fieldnames) # use the DictWriter module and declare your fieldnames
                    f_out.writeheader() # write the header
                    for elem in genes: # iterate over your structure (genes) and for every element write a row
                        elem.window_size = window # update the attribute of window sie to be the current window just before you write
                        f_out.writerow(elem.__dict__) # elem.__dict__ is how we access all of the attributes of an instance
        window += increment




        reset(genes)
#===============================================
# Reset the Deque for the next window
#===============================================
def reset(deque):
    """Resets the attributes altered during the course of the algorithm for use in the next window."""
    for elem in deque:

        elem.DNA_left = 0
        elem.DNA_intra = 0
        elem.DNA_right = 0

        elem.LTR_left = 0
        elem.LTR_intra = 0
        elem.LTR_right = 0

        elem.LINE_left = 0
        elem.LINE_intra = 0
        elem.LINE_right  = 0

        elem.Unknown_left = 0
        elem.Unknown_intra = 0
        elem.Unknown_right = 0

        elem.MULE_left = 0
        elem.MULE_intra = 0
        elem.MULE_right = 0

        elem.Gypsy_left = 0
        elem.Gypsy_intra = 0
        elem.Gypsy_right = 0

        elem.Unknown_fam_left = 0
        elem.Unknown_fam_intra = 0
        elem.Unknown_fam_right = 0

        elem.CMC_EnSpm_left = 0
        elem.CMC_EnSpm_intra = 0
        elem.CMC_EnSpm_right = 0

        elem.Copia_left = 0
        elem.Copia_intra = 0
        elem.Copia_right = 0

        elem.LINE_fam_left = 0
        elem.LINE_fam_intra = 0
        elem.LINE_fam_right = 0

        elem.hAT_left = 0
        elem.hAT_intra = 0
        elem.hAT_right  = 0

        elem.PIF_Harbinger_left  = 0
        elem.PIF_Harbinger_intra  = 0
        elem.PIF_Harbinger_right  = 0

        elem.LTRUnknown_left = 0
        elem.LTRUnknown_intra = 0
        elem.LTRUnknown_right = 0

        elem.they_are_inside = 0

#===================================================================================
# Here are my algorithm's sorting functions
#===================================================================================
def left_outside_from_window(window,left_window_stop,left_window_start,left_density):
    """TE comes from outside the window, and ends within the gene, ths part deals with the window to the gene start. Similar to right_outside_from_window."""
    outside_chunk = left_window_stop - left_window_start  + 1
    left_density += outside_chunk/window
    return left_density

def left_inside(t_stop,g_start,g_stop,t_type,intra_density):
    """Coming from the left, this deals the part of the TE that is inside the gene."""
    inside_chunk = t_stop - g_start + 1
    intra_density += inside_chunk/(g_stop - g_start + 1) # divide the inside part by g_length
    return intra_density

def left_outside_in_window(window,left_window_stop,t_start,left_density):
    """ TE starts within the window but ends inside the gene. Only deals with the portion that is in the window. """
    outside_chunk = left_window_stop - t_start + 1
    left_density += outside_chunk/window
    return left_density

def right_inside(t_start,g_start,g_stop,t_type,intra_density):
    """right_inside() calculates the intra-gene density for a TE."""
    inside_chunk = g_stop - t_start + 1
    intra_density += inside_chunk/(g_stop - g_start + 1)
    return intra_density

def right_outside_in_window(window,right_window_start,t_stop,right_density):
    """ Calculates outside of gene density when the TE ends within the window."""
    outside_chunk = t_stop - right_window_start + 1
    right_density += outside_chunk/window
    return right_density

def right_outside_from_window(window,right_window_start,right_window_stop,right_density):
    """ Outside the gene, extends past the window, the entire inside of the window is covered """
    outside_chunk = right_window_stop - right_window_start + 1
    right_density += outside_chunk/window
    return right_density

def L_prox_check(left_window_stop,t_stop,prox_left):
    t_stop += 1
    distance = left_window_stop - t_stop
    if prox_left == None:
        prox_left = distance
    if distance < prox_left:
        prox_left = distance
    return prox_left

def R_prox_check(right_window_start,t_start,prox_right):
    t_start -= 1
    distance = t_start - right_window_start
    if prox_right == None:
        prox_right = distance
    if distance < prox_right:
        prox_right = distance
    return prox_right

#===================================================================================
if __name__ == '__main__':
    pass
