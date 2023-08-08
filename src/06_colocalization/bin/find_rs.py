import sys
import os
path_rsid=sys.argv[1]
path_3qtl=sys.argv[2]
files_rsid = os.listdir(path_rsid)
files_3qtl = os.listdir(path_3qtl)
dict_rsid1={}
dict_rsid2={}
dict_rsid3={}
dict_rsid4={}
dict_rsid5={}
dict_rsid6={}
dict_rsid7={}
dict_rsid8={}
dict_rsid9={}
dict_rsid10={}
dict_rsid11={}
dict_rsid12={}
dict_rsid13={}
dict_rsid14={}
dict_rsid15={}
dict_rsid16={}
dict_rsid17={}
dict_rsid18={}
dict_rsid19={}
dict_rsid20={}
dict_rsid21={}
dict_rsid22={}
for file in files_rsid:
	if not os.path.isdir(file):
		f = open(path_rsid+"/"+file,"r")
		file_split=file.split('_')
		chrnum=file_split[0]
		iter_f=iter(f)
		if chrnum == "chr1":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid1[key]=value
		if chrnum == "chr2":
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid2[key]=value
		if chrnum == "chr3":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid3[key]=value
		if chrnum == "chr4":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid4[key]=value
		if chrnum == "chr5":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid5[key]=value

		if chrnum == "chr6":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid6[key]=value
		if chrnum == "chr7":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid7[key]=value
		if chrnum == "chr8":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid8[key]=value
		if chrnum == "chr9":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid9[key]=value
		if chrnum == "chr10":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid10[key]=value
		if chrnum == "chr11":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid11[key]=value
		if chrnum == "chr12":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid12[key]=value
		if chrnum == "chr13":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid13[key]=value
		if chrnum == "chr14":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid14[key]=value
		if chrnum == "chr15":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid15[key]=value
		if chrnum == "chr16":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid16[key]=value
		if chrnum == "chr17":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid17[key]=value
		if chrnum == "chr18":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid18[key]=value
		if chrnum == "chr19":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid19[key]=value
		if chrnum == "chr20":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid20[key]=value
		if chrnum == "chr21":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid21[key]=value
		if chrnum == "chr22":	
			for line in iter_f:
				line_split=filter(None, line.split('\t'))
				if line_split[0].isdigit():
					key=line_split[1]
					value=line_split[3]
					dict_rsid22[key]=value
		f.close()
		print("The dictionary of "+file+" is done")			
		print("..................................")
for file in files_3qtl:
	if not os.path.isdir(file):
		f = open(path_3qtl+"/"+file,"r")
		file_split = file.split('.')
		tissue=file_split[0]
		w = open(path_3qtl+"/../3QTLs_rsid/"+tissue+"_rsid.txt","w+")
		iter_f=iter(f)
		for line in iter_f:
			line_split=filter(None, line.split('\t'))
			if line_split[0] != "chr":
				key=line_split[1].split("_")[1]
				if line_split[0] == "chr1" and key in dict_rsid1:
					line=line[:-1]+"\t"+dict_rsid1[key]+"\n"
					w.write(line)
				if line_split[0] == "chr2" and key in dict_rsid2:
					line=line[:-1]+"\t"+dict_rsid2[key]+"\n"
					w.write(line)
				if line_split[0] == "chr3" and key in dict_rsid3:
					line=line[:-1]+"\t"+dict_rsid3[key]+"\n"
					w.write(line)
				if line_split[0] == "chr4" and key in dict_rsid4:
					line=line[:-1]+"\t"+dict_rsid4[key]+"\n"
					w.write(line)
				if line_split[0] == "chr5" and key in dict_rsid5:
					line=line[:-1]+"\t"+dict_rsid5[key]+"\n"
					w.write(line)
				if line_split[0] == "chr6" and key in dict_rsid6:
					line=line[:-1]+"\t"+dict_rsid6[key]+"\n"
					w.write(line)
				if line_split[0] == "chr7" and key in dict_rsid7:
					line=line[:-1]+"\t"+dict_rsid7[key]+"\n"
					w.write(line)
				if line_split[0] == "chr8" and key in dict_rsid8:
					line=line[:-1]+"\t"+dict_rsid8[key]+"\n"
					w.write(line)
				if line_split[0] == "chr9" and key in dict_rsid9:
					line=line[:-1]+"\t"+dict_rsid9[key]+"\n"
					w.write(line)
				if line_split[0] == "chr10" and key in dict_rsid10:
					line=line[:-1]+"\t"+dict_rsid10[key]+"\n"
					w.write(line)
				if line_split[0] == "chr11" and key in dict_rsid11:
					line=line[:-1]+"\t"+dict_rsid11[key]+"\n"
					w.write(line)
				if line_split[0] == "chr12" and key in dict_rsid12:
					line=line[:-1]+"\t"+dict_rsid12[key]+"\n"
					w.write(line)
				if line_split[0] == "chr13" and key in dict_rsid13:
					line=line[:-1]+"\t"+dict_rsid13[key]+"\n"
					w.write(line)
				if line_split[0] == "chr14" and key in dict_rsid14:
					line=line[:-1]+"\t"+dict_rsid14[key]+"\n"
					w.write(line)
				if line_split[0] == "chr15" and key in dict_rsid15:
					line=line[:-1]+"\t"+dict_rsid15[key]+"\n"
					w.write(line)
				if line_split[0] == "chr16" and key in dict_rsid16:
					line=line[:-1]+"\t"+dict_rsid16[key]+"\n"
					w.write(line)
				if line_split[0] == "chr17" and key in dict_rsid17:
					line=line[:-1]+"\t"+dict_rsid17[key]+"\n"
					w.write(line)
				if line_split[0] == "chr18" and key in dict_rsid18:
					line=line[:-1]+"\t"+dict_rsid18[key]+"\n"
					w.write(line)
				if line_split[0] == "chr19" and key in dict_rsid19:
					line=line[:-1]+"\t"+dict_rsid19[key]+"\n"
					w.write(line)
				if line_split[0] == "chr20" and key in dict_rsid20:
					line=line[:-1]+"\t"+dict_rsid20[key]+"\n"
					w.write(line)
				if line_split[0] == "chr21" and key in dict_rsid21:
					line=line[:-1]+"\t"+dict_rsid21[key]+"\n"
					w.write(line)
				if line_split[0] == "chr22" and key in dict_rsid22:
					line=line[:-1]+"\t"+dict_rsid22[key]+"\n"
					w.write(line)
			else:
				line=line[:-1]+"\t"+"rsid"+"\n"
				w.write(line)
		w.close()			
		f.close()
		print("All the rsid in "+file+" is found")
		print(".................................")

					
		
