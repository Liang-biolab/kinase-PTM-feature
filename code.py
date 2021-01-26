import pickle
import os
import sys


dd = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 
	 'GLY': 'G', 'GLN': 'Q', 'GLU': 'E', 'HIS': 'H', 'ILE': 'I',
	 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PRO': 'P', 'PHE': 'F',
	 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'TRP': 'W', 'VAL': 'V'}
standard_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLY', 'GLN', 'GLU', 'HIS', 'ILE', 
						'LEU', 'LYS', 'MET', 'PRO', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


def get_pdb_sequence(pdb_file, chain):
	# 从pdb文件中得到指定的链的序列
	pdb_sequence  = ''
	with open(pdb_file, 'r') as f:
		for line in f:
			if line[:6] != 'SEQRES':
				continue
			if line[11] != chain:
				continue
			line = line.strip().split(' ')
			while ' ' in line:
				line.remove(' ')
			while '' in line:
				line.remove('')
			for i in range(4, len(line)):
				aa = line[i]
				if aa not in standard_amino_acids:
					continue
				pdb_sequence += dd[aa]

	return pdb_sequence


def get_uniprot_sequence(uniprot_id):
	uniprot_sequence = ''
	with open('uniprot.fasta', 'r') as f:
		while  True:
			line = f.readline()
			if line == '':  # 如果到了文件末尾 则退出循环
				break
			# print(line)
			if line[0] == '>':
				uniprotid = line[1:].strip()
			elif uniprotid == uniprot_id:   # 如果该swissid 属于我要找的，就保存它的seq
				line = line.strip()
				uniprot_sequence += line

	return uniprot_sequence


def write_fasta(uniprot_sequence, pdb_sequence, uniprot, pdb):
	with open('./example/' + uniprot + '_' + pdb + '.fasta', 'w') as f:
		f.write('>' + uniprot + '\n')
		f.write(uniprot_sequence + '\n')

		f.write('>' + pdb + '\n')
		f.write(pdb_sequence)

	
def run_clustal(uniprot, pdb):
	input_file = './example/' + uniprot + '_' + pdb + '.fasta'
	output_file = './example/' + uniprot + '_' + pdb + '.fasta.aln'  # clusta的结果文件

	cmdr = os.system('clustalo.exe -i %s --seqtype=Protein -o %s'% (input_file, output_file))

	os.remove(input_file)


def get_pdb_msa_seq(uniprot, pdb):
	'''
	读取pdb序列和swiss序列序列比对后的结果文件，提取出pdb的结果序列
	'''

	file = './example/' + uniprot + '_' + pdb + '.fasta.aln' 
	pdb_seq = ''
	
	with open(file, 'r') as f:
		find_pdb = False
		for line in f:
			line = line.strip()
			if line[0] == '>':
				if line.find(pdb) != -1:
					find_pdb = True
					continue
			if find_pdb == True:
				pdb_seq += line

	os.remove(file)
				
	return pdb_seq


def get_pdb_interval(pdb_seq):
	s = pdb_seq
	p = pdb_seq.split('-')
	while '' in p:
		p.remove('')
	temp_interval = []
	left = 1
	for pp in p:
		p_in_s = s.find(pp) + 1
		temp_interval.append((left, left+len(pp)-1, p_in_s, p_in_s+len(pp)-1))
		left += len(pp)

	with open('./example/pdb_swiss.txt', 'w') as f:
		f.write('pdb' + '\t' + 'uniprot' + '\n')
		for interval in temp_interval:
			start = interval[0]
			end = interval[1]
			count = interval[2]
			for i in range(start, end+1):
				f.write(str(i) + '\t' + str(count) + '\n')
				count += 1
		

		# break

		
if __name__ == '__main__':
	pdb = '1a4m'
	pdb_file = './example/' + pdb + '.pdb'
	chain = 'A'
	uniprot = 'P03958'

	pdb_sequence = get_pdb_sequence(pdb_file, chain)
	uniprot_sequence = get_uniprot_sequence(uniprot)

	# 将uniprot_sequence 和 pdb_sequence写入fasta，以便接下来进行多序列比对
	# 保存的文件名为 uniprot_pdb.fasta
	write_fasta(uniprot_sequence, pdb_sequence, uniprot, pdb)

	# 使用clustal进行序列比对，结果文件名为 uniprot_pdb.fasta.aln
	run_clustal(uniprot, pdb)

	pdb_sequence_new = get_pdb_msa_seq(uniprot, pdb)

	get_pdb_interval(pdb_sequence_new) # 最后的结果文件为pdb_uniprot.txt
	
