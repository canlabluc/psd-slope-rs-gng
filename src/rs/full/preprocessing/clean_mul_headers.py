import os
import glob

## Parameters
src_model = '' # dmn, frontal, dorsal, ventral
import_path = ''
export_path = ''


def get_filelist(import_path):
    matfiles = []
    for root, dirs, files in os.walk(import_path):
        matfiles += glob.glob(os.path.join(root, '*.mul'))
    return matfiles

file_list = get_filelist(import_path)
for file in file_list:
    newfile_name = export_path + '/' + file.split('/')[-1]
    oldfile = open(file, 'r')
    newfile = open(newfile_name, 'w')
    lines = oldfile.readlines()
    if src_model == 'dmn':
        lines[1] = ' PCC PCCr PCCv PCCh mPFC mPFCr mPFCv mPFCh LAG LAGr LAGv LAGh RAG RAGr RAGv RAGh LLatT LLatTe1 LLatTe2 LLatTe3 RLatT RLatTe1 RLatTe2 RLatTe3 Noise1L Noise1L1 Noise1L2 Noise1L3 Noise1R Noise1R1 Noise1R2 Noise1R3 Noise2L Noise2L1 Noise2L2 Noise2L3 Noise2R Noise2R1 Noise2R2 Noise2R3 Noise1M Noise1M1 Noise1M2 Noise1M3 Noise2M Noise2M1 Noise2M2 Noise2M3'
    elif src_model == 'frontal':
        lines[1] = ' LdlPFC LdlPFC1 LdlPFC2 LdlPFC3 RdlPFC RdlPFC1 RdlPFC2 RdlPFC3 LFRont LFront1 LFront2 LFront3 RFront RFront1 RFront2 RFront3 LIPL LIPLr LIPLv LIPLh RIPL RIPLr RIPLv RIPLh LIPS LIPSr LIPSv LIPSh RIPS RIPSr RIPSv RIPSh Noise1L Noise1L1 Noise1L2 Noise1L3 Noise1R Noise1R1 Noise1R2 Noise1R3 Noise2L Noise2L1 Noise2L2 Noise2L3 Noise2R Noise2R1 Noise2R2 Noise2R3 NoiseF NoiseF1 NoiseF2 NoiseF3'
    elif src_model == 'dorsal':
        lines[1] = ' LFEF LFEFr LFEFv LFEFh RFEF RFEFr RFEFv RFEFh LaIPS LaIPSr LaIPSv LaIPSh RaIPS RaIPSr RaIPSv RaIPSh LpIPS LpIPSr LpIPSv LpIPSh RpIPS RpIPSr RpIPSv RpIPSh Noise1L Noise1L1 Noise1L2 Noise1L3 Noise1R Noise1R1 Noise1R2 Noise1R3 Noise2L Noise2L1 Noise2L2 Noise2L3 Noise2R Noise2R1 Noise2R2 Noise2R3 Noise3L Noise3L1 Noise3L2 Noise3L3 Noise4R Noise4R1 Noise4R2 Noise4R3'
    elif src_model == 'ventral':
        lines[1] = ' LIFG LIFGr LIFGv LIFGh RIFG RIFGr RIFGv RIFGh LMFG LMFGr LMFGv LMFGh RMFG RMFGr RMFGv RMFGh LTPJ LTPJr LTPJv LTPJh RTPJ RTPJr RTPJv RTPJh LSTG LSTGr LSTGv LSTGh RSTG RSTGr RSTGv RSTGh NoiseL NoiseL1 NoiseL2 NoiseL3 NoiseR NoiseR1 NoiseR2 NoiseR3 NoiseF NoiseF1 NoiseF2 NoiseF3 Noise Noise1 Noise2 Noise3'

    newfile.writelines(lines)

    print("Finished: ", newfile_name)

    oldfile.close()
    newfile.close()


