import sys
import os
from subprocess import *
from multiprocessing import Pool
from optparse import OptionParser
from Bio import SeqIO

#
BWA_PATH="bwa"
SAMTOOLS_PATH="samtools"
WORKING_FOlDER="./out/"
CLIP_CUTOFF=13
S_DELIM="~"
BRKPNT_EXTND=1000


#load in the regions from file to a list
def load_regions(sf_region):#
    global S_DELIM
    l_regions_db = []
    with open(sf_region) as fin_region:
        #cnt=1
        for line in fin_region:
            fields=line.split()
            chrm=fields[0]
            start=int(fields[1])
            end=int(fields[2])
            l_regions_db.append((chrm, start, end))
    return l_regions_db


##parse the fasta sequences from reference
def gnrt_region_ref(sf_ref, l_regions_db, i_extd):
    global WORKING_FOlDER
    for record in SeqIO.parse(sf_ref, "fasta"):
        cur_chr=str(record.id)
        for reg in l_regions_db:
            chr=reg[0]
            start=reg[1] ##extend the breakpoint for some specific distance
            end=reg[2]
            rgn_start=start-i_extd
            rgn_end=end+i_extd

            if chr==cur_chr:
                aseq=str(record.seq[rgn_start:rgn_end+1])
                sid = chr + S_DELIM + str(start) + S_DELIM + str(end)
                new_file = WORKING_FOlDER + sid + ".fa"
                with open(new_file,"w") as fout_new:
                    fout_new.write(">"+sid+"\n")
                    fout_new.write(aseq+"\n")

def run_cmd(cmd):
    Popen(cmd, shell = True, stdout = PIPE).communicate()


##index the regions fa files in parallel
def index_region_ref(m_regions_db, n_jobs):
    global WORKING_FOlDER, BWA_PATH

    cmd_list=[]
    for reg in m_regions_db:
        chr = reg[0]
        start = reg[1]  ##extend the breakpoint for some specific distance
        end = reg[2]
        sid = chr + S_DELIM + str(start) + S_DELIM + str(end)
        sf_region=WORKING_FOlDER + sid + ".fa"

        if os.path.exists(sf_region)==False:
            print "File {0} does not exist!!".format(sf_region)
            return
        cmd="{0} index {1}".format(BWA_PATH, sf_region)
        cmd_list.append(cmd)

    pool = Pool(n_jobs)
    pool.map(run_cmd, cmd_list, 1)
    pool.close()
    pool.join()


##Check whether reads is qualified clipped reads
def is_qualified_clipped(cigar, cutoff_len):
    l=len(cigar)
    signal=[]
    lenth=[]
    temp=""
    for i in range(l):
        if cigar[i]>="0" and cigar[i]<="9":
            temp=temp+cigar[i]
        else:
            signal.append(cigar[i])
            try:
                lenth.append(int(temp))
            except ValueError:
                print "Error: ", cigar, temp
            temp=""

    b_hardclip=False
    cnt_m=0
    clip_flag=0
    left_clip_len=0
    right_clip_len=0

    for i in range(len(signal)):
        if signal[i]=="M":
            cnt_m=cnt_m+lenth[i]

    if (signal[0]=="S" or signal[0]=="H") and lenth[0]>=cutoff_len:#left-clip
        clip_flag=1
        left_clip_len=lenth[0]
        if signal[0]=="H":
            b_hardclip=True
    if (signal[len(signal)-1]=="S" or signal[len(signal)-1]=="H") and lenth[len(signal)-1]>=cutoff_len: #right-clip
        clip_flag=clip_flag+2 #if this is 3, then both side clipped
        right_clip_len=lenth[len(signal)-1]
        if signal[len(signal)-1]=="H":
            b_hardclip=True

    return b_hardclip, cnt_m, clip_flag, left_clip_len, right_clip_len


def is_fully_map(cigar, max_clip_len):
    l=len(cigar)
    signal=[]
    lenth=[]
    temp=""
    for i in range(l):
        if cigar[i]>="0" and cigar[i]<="9":
            temp=temp+cigar[i]
        else:
            signal.append(cigar[i])
            try:
                lenth.append(int(temp))
            except ValueError:
                print "Error: ", cigar, temp
            temp=""

    cnt_m=0
    total_len=0
    for i in range(len(signal)):
        if signal[i]=="M":
            cnt_m=cnt_m+lenth[i]
        total_len=total_len+lenth[i]

    if (total_len-cnt_m)<=max_clip_len:
        return True, cnt_m
    return False, 0


#check whether a pair is discordant or not
#TBD: Note, more strict conditions: only keep the discordant pairs related to the event type
def is_discordant(flag, mean_is, std_is, tlen):
    # if flag&2!=0:#0x2 is set, then each segments are properly mapped
    #     return False
    # else:
    #     return True
    if flag&8!=0:#read mapped and mate unmapped, then not an discordant pair
        return False

    if float(tlen) > (mean_is+3.0*std_is) or float(tlen) < (mean_is-3.0*std_is):
        return True
    else:
        return False



#parse the alignment in the given region
#output:
#1. the clipped parts in fastq format
#2. the clip / discordant reads in dictionary
def parse_alignment_gnrt_disc_clip(sf_sam, m_insert_size, sf_clip):
    global CLIP_CUTOFF
    m_clip={} #[qname, ((rname, map_pos, cigar), flag, (clip_pos ,rnext, pnext))]
    m_discord={} ##each record in format [qname, ((rname, map_pos, cigar), flag, (tlen ,rnext, pnext))]
    with open(sf_sam) as fin_sam, open(sf_clip, "w") as fout_clip:
        cnt=0
        for line in fin_sam:
            if line[0]=="@":
                continue
            cnt=cnt+1
            fields=line.split()
            qname=fields[0]
            flag=int(fields[1])
            rname=fields[2]
            cigar=fields[5]
            map_pos=int(fields[3])#mapping position
            rnext=fields[6]
            pnext=int(fields[7])
            tlen=abs(int(fields[8]))

            if cigar=="*":#unmapped are not interesting
                continue

            if rnext != "=":
                m_discord[qname] = ((rname, map_pos, cigar), flag, (tlen, rnext, pnext))
            else:
                ##check discordant reads
                #get the read group information from the alignment
                rg_id=""
                for info in fields[11:]:
                    if len(info)>3 and info[:3]=="RG:":
                        rg_id=info[5:]
                if rg_id in m_insert_size:
                    mean_is=float(m_insert_size[rg_id][0])
                    std_is=float(m_insert_size[rg_id][1])
                    if is_discordant(flag, mean_is, std_is, tlen)==True:
                        m_discord[qname] = ((rname, map_pos, cigar), flag, (tlen, rnext, pnext))
                else:
                    print "Read group ", rg_id, " not found in the library!!!!!"

            #check clipped reads
            b_hardclip, cnt_m, clip_flag, left_clip_len, right_clip_len=is_qualified_clipped(cigar, CLIP_CUTOFF)
            if b_hardclip==True or clip_flag==0:#hard-clip or fully mapped reads are not interested
                continue
            if left_clip_len>0 and right_clip_len>0:#both end clipped reads are not interested
                continue

            read_seq=fields[9]
            quality_seq=fields[10]
            clipped_seq=""
            clipped_qseq=""
            clip_pos=map_pos
            bleft_clip=1
            read_len=len(read_seq)
            if right_clip_len>0:#right clip
                clip_pos=map_pos+cnt_m-1
                clipped_seq=read_seq[read_len-right_clip_len:]
                clipped_qseq=quality_seq[read_len-right_clip_len:]
                bleft_clip=0
            else:#left clip
                clipped_seq=read_seq[:left_clip_len]
                clipped_qseq=quality_seq[:left_clip_len]

            # if clip_pos not in m_clip_pos:
            #     m_clip_pos[clip_pos]=0
            # m_clip_pos[clip_pos]=m_clip_pos[clip_pos]+1

            m_clip[qname]=((rname, map_pos, cigar), flag, (clip_pos, rnext, pnext))

            ##parse out the clipped part and save into sf_clip
            fout_clip.write("@"+qname+"_"+str(cnt)+"_"+str(bleft_clip)+"\n")
            fout_clip.write(clipped_seq+"\n")
            fout_clip.write("+\n")
            fout_clip.write(clipped_qseq+"\n")
    return m_discord, m_clip



##this function collect the clip positions from the alignments
##Note that, all the pos here are "absolute positions" on this new ref, so need to plus the "start position" in the old "ref"
##Also, need to call the "clip position" according to the clip direction (left or right) of the original clip reads
##  this information saved in the last char of the read id (0 indicates original read is right clip, 1 indicate left clip)
## Here give a "hard code" clip length 4, if smaller than this, then still considered as fully mapped
def parse_clipped_algmt(sf_sam, start_pos):
    m_clip={}
    max_clip_len=4 ##################################################################################
    with open(sf_sam) as fin_sam:
        for line in fin_sam:
            if line[0]=="@":
                continue
            fields=line.split()
            qname=fields[0]
            flag=int(fields[1])
            cigar=fields[5]
            map_pos=int(fields[3])#mapping position
            if cigar=="*":#unmapped are not interesting
                continue

            bfull_map, cnt_m=is_fully_map(cigar, max_clip_len)
            if bfull_map==False:
                continue

            bleft_clip=True
            if qname[-1]=="0":
                bleft_clip=False

            clip_pos=map_pos + start_pos -1 #####################################################here neeed to -1?????????????????????????????????
            if bleft_clip==True:
                clip_pos=clip_pos + cnt_m

            qfields=qname.split("_")
            qqname=qfields[0]
            m_clip[qqname]=(flag, clip_pos)

            # if clip_pos not in m_clip_pos:
            #     m_clip_pos[clip_pos]=0
            # m_clip_pos[clip_pos]=m_clip_pos[clip_pos]+1
    return m_clip



#mpart1 in format: [key, ((rname, map_pos, cigar), flag, (clip_pos, rnext, pnext))]
#mpart2 in format: [key, (flag, clip_pos)]
# check two consistents:
# 1) clip position are within specific distance from the breakpoints (by checking against "dist_clip_pos")
# 2) two parts of the read should have consistent orientation (by checking "b_expect_ori")
##
def check_clip_consistent(mpart1, mpart2, start, end, b_expect_ori, dist_clip_pos):
    m_clip={}
    for qname in mpart1:
        if mpart2.has_key(qname)==False:
            continue

        flag1=int(mpart1[qname][1])
        clip_pos1=int(mpart1[qname][2][0])#
        b_rc1=False
        if flag1 & 16 != 0:
            b_rc1=True

        flag2=int(mpart2[qname][0])#
        clip_pos2=int(mpart2[qname][1])#
        b_rc2=False
        if flag2 & 16 != 0:
            b_rc2=True

        #
        b_l_d = abs(start-clip_pos1) < dist_clip_pos
        b_r_d = abs(end-clip_pos2) < dist_clip_pos
        b_within_dist1= (b_l_d and b_r_d)


        b_l_d2 = abs(end-clip_pos1) < dist_clip_pos
        b_r_d2 = abs(start-clip_pos2) < dist_clip_pos
        b_within_dist2= (b_l_d2 and b_r_d2)

        b_ori = (b_rc1==b_rc2)

        if (b_ori==b_expect_ori) and (b_within_dist1 or b_within_dist2): #same direction, and "clip_pos" is within specific distance from the "breakpoints"
            m_clip[qname]=b_rc1

    return m_clip

#
#
#

##Note: set "hard code" window size 10
def run_parse_region(regions):
    global SAMTOOLS_PATH
    global BWA_PATH
    global WORKING_FOlDER
    global CLIP_CUTOFF
    global BRKPNT_EXTND

    #window_size=10 ############################################hard code here #########################################
    dist_clip_pos=50 ###################################################################################################
    b_expect_ori=True ##same direction #######################################
    reg=regions[0] #(s,i,i)
    sf_bam=regions[1]
    sf_insert_size=regions[2]
    m_insert_size = load_insert_size(sf_insert_size)  ##insert size dictionary


    chr_a=reg[0]
    start_a=reg[1]
    start_pos = start_a - BRKPNT_EXTND
    end_a=reg[2]
    end_pos=end_a+BRKPNT_EXTND

    sreg=chr_a+S_DELIM+str(start_a)+S_DELIM+str(end_a)
    sf_sam=WORKING_FOlDER+sreg+".sam"
    cmd="{0} view {1} {2}:{3}-{4} > {5}".format(SAMTOOLS_PATH, sf_bam, chr_a, start_pos, end_pos, sf_sam)
    run_cmd(cmd)

    sf_clip=WORKING_FOlDER+sreg+".clip"
    m_discord, m_clip=parse_alignment_gnrt_disc_clip(sf_sam, m_insert_size, sf_clip)

    #align the clipped part to specific region
    sref=WORKING_FOlDER+sreg+".fa"
    sf_clip_align= WORKING_FOlDER + sreg+"_align"+".sam"
    cmd="{0} mem -T {1} {2} {3} > {4}".format(BWA_PATH, CLIP_CUTOFF, sref, sf_clip, sf_clip_align)
    run_cmd(cmd)

    #start_pos=start_a-BRKPNT_EXTND
    m_clip_part_algmt=parse_clipped_algmt(sf_clip_align, start_pos)

    ##check the two parts of the clipped reads
    m_clip_consistent=check_clip_consistent(m_clip, m_clip_part_algmt, start_a, end_a, b_expect_ori, dist_clip_pos)

    #output the discordant and clipped reads for the event
    sf_rslt=WORKING_FOlDER+"{0}.rslt".format(sreg)
    with open(sf_rslt,"w") as fout_rslt:
        for qname in m_discord:
            fout_rslt.write(qname+"\tdisc\n")
        for qname in m_clip_consistent:
            fout_rslt.write(qname+"\tclip\n")



##ToDoList: 1. check the orientation between discordant pairs and clippe parts whether are consistent or not!!!!
##ToDoList: 2. check the orientation of the discordant pairs against the SV type
##TodoList: 3. merge the bam

#parse in parallel mode
def parse_disc_clip_reads(sf_algnmt, l_regions, sf_insert_size, n_jobs):
    global WORKING_FOlDER
    reg_list=[]
    for reg in l_regions:# reg is a tuple: (s,i,i)
        region=(reg, sf_algnmt, sf_insert_size)
        reg_list.append(region)

    pool = Pool(n_jobs)
    pool.map(run_parse_region, reg_list, 1)
    pool.close()
    pool.join()

    #merge the results
    #cmd="cat {0}*.rslt > {1}disc_clip.reads".format(WORKING_FOlDER, WORKING_FOlDER)
    #run_cmd(cmd)

#remove the temporary files

def load_insert_size(sf_insert_size):
    m_insert_size={}
    with open(sf_insert_size) as fin_insert_size:
        for line in fin_insert_size:
            fields=line.split()
            rg_name=fields[0]
            mean_insert_size=float(fields[1])
            std_insert_size=float(fields[2])
            m_insert_size[rg_name]=(mean_insert_size, std_insert_size)
    return m_insert_size


#Retrieve the discordant and clipped reads from alignment
def retrieve_reads(sf_ref, sf_algnmt, sf_variation, sf_insert_size, iextd, n_jobs, sf_out_prefix):
    global WORKING_FOlDER
    global BRKPNT_EXTND
    WORKING_FOlDER=sf_out_prefix
    BRKPNT_EXTND=iextd

    if os.path.exists(WORKING_FOlDER)==False:
        print "Working folder", WORKING_FOlDER, "doesn't exist!!!"
        return

    if sf_out_prefix[-1]!="/":
        WORKING_FOlDER = sf_out_prefix+"/"

    l_regions_db=load_regions(sf_variation) ####[(s,i,i)] ###
    if len(l_regions_db)<1:
        print "The variation file is empty!!"
        return
    gnrt_region_ref(sf_ref, l_regions_db, iextd)
    index_region_ref(l_regions_db, n_jobs)

    m_insert_size=load_insert_size(sf_insert_size)##insert size dictionary
    if len(m_insert_size)<1:
        print "Please provide the insert size information for each group!!!"
        return
    parse_disc_clip_reads(sf_algnmt, l_regions_db, sf_insert_size, n_jobs)



##parse the options
def parse_option():
    parser = OptionParser()
    # parser.add_option("-s", "--select",
    #                   action="store_true", dest="select", default=False,
    #                   help="run select: select subcontigs from give contig file")
    parser.add_option("-t", "--retrive",
                      action="store_true", dest="retrieve", default=False,
                      help="run retrive pipeline")
    parser.add_option("-r", "--ref", dest="reference",
                      help="Reference file that will be aligned to", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="The original bam file (sorted and indexed) ", metavar="FILE")
    parser.add_option("-v", "--var", dest="variation",
                      help="The variation file in bed format (first 3 columns are required) ", metavar="FILE")
    parser.add_option("-i", "--is", dest="insert_size",
                      help="Insert size information) ", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-e", "--lenth", dest="extlen", type="int",
                      help="breakpoint extend length")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of threads")

    (options, args) = parser.parse_args()
    return (options, args)



##main function
if __name__ == '__main__':

    (options, args) = parse_option()

    if options.retrieve:
        sf_ref=options.reference
        sf_variation=options.variation
        sf_algnmt=options.bam
        i_extend=options.extlen
        n_jobs=options.cores
        sf_out_prefix=options.output
        sf_insert_size=options.insert_size

        #retrieve the interested reads
        retrieve_reads(sf_ref, sf_algnmt, sf_variation, sf_insert_size, i_extend, n_jobs, sf_out_prefix)

    else:
        print "Wrong parameters!"


####