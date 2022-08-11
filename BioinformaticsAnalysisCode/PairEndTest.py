import gzip

fq1_path = "C:\\Users\\Jeff\\OneDrive\\桌面\\BioinformaticsAnalysisCode\\W0gDNA.gw8around_R1.fastq.gz"
f = gzip.open(fq1_path, 'rb')
fq1_name=[]
for line in f.readlines(): # 按行进行读取
    s = line.decode() # 读取之后要进行解码
    if s.startswith("@FP"):
        fq1_name.append(s.split("/1")[0]) # s 为string类型，就是我们读取的文件中的一行

with open("C:\\Users\\Jeff\\OneDrive\\桌面\\BioinformaticsAnalysisCode\\fq1.txt", 'w') as f:
    for i in fq1_name:
        f.write(i+"\n")



fq2_path = "C:\\Users\\Jeff\\OneDrive\\桌面\\BioinformaticsAnalysisCode\\W0gDNA.gw8around_R2.fastq.gz"
g = gzip.open(fq2_path, 'rb')
fq2_name=[]
for line in g.readlines(): # 按行进行读取
    s = line.decode() # 读取之后要进行解码
    if s.startswith("@FP"):
        fq2_name.append(s.split("/2")[0]) # s 为string类型，就是我们读取的文件中的一行
with open("C:\\Users\\Jeff\\OneDrive\\桌面\\BioinformaticsAnalysisCode\\fq2.txt", 'w') as g:
    for i in fq2_name:
        g.write(i+"\n")


with open("C:\\Users\\Jeff\\OneDrive\\桌面\\BioinformaticsAnalysisCode\\intersect.txt", 'w') as g:
    for i in list(set(fq1_name) & set(fq2_name)):
        g.write(i+"\n")
