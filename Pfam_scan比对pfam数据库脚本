Pfam_scan比对pfam数据库
#############################################################
####################pfam_scan 注释###########################
######################LorMe 实验室############################
######################    zyz    ############################
#############################################################
一、软件说明
################################################
# pfam数据库简介:
# Pfam是一个蛋白家族数据库，其中Pfam-A是手工确定的高质量的蛋白家族，Pfam-B是自动注释的，是对A的补充。目前已更新到34.0,包含超过13000个手工确定的蛋白家族。任选一版本即可，需要两个文件，Pfam-A.hmm.gz和 Pfam-A.hmm.dat.gz。
# 其中两个数据库:
1.高质量，手工确定的Pfam-A，
2.自动注释的Pfam-B数据库(该部分数据产生是根据ADDA算法，是对A的补充)
# 一般Pfam-B数据库中数据质量不高，因此主要选择Pfam-A数据库对蛋白序列进行注释。基于蛋白序列获取蛋白domain的原理是参考的blast原理，当蛋白序列相似度高时，其domain也相似，因此，相当于是基于序列预测结构。
# 但这个结构毕竟是预测的，对其domain的可靠性要取cutoff。一般情况下为了保证正确率，e-value值的cutoff定为10-6，(凌恩公司设置为10-5)这一般视研究课题而定。
################################################
# 参考网站：
# PfamScan及fam数据库安装:http://blog.sina.com.cn/s/blog_6e18b79b0102wyvl.html
# 数据库Pfam在线注释以及本地化全攻略:https://www.yunbios.net/Pfam.html
# Pfam官网:http://pfam.xfam.org/
# 数据库下载地址:http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
在线注释网站:https://www.ebi.ac.uk/Tools/pfa/pfamscan/
################################################
# 必读说明：
# 利用PfamScan寻找同源基因家族:https://blog.csdn.net/weixin_43840576/article/details/86589573?utm_medium=distribute.pc_relevant.none-task-blog-2%7Edefault%7EBlogCommendFromMachineLearnPai2%7Edefault-4.control&depth_1-utm_source=distribute.pc_relevant.none-task-blog-2%7Edefault%7EBlogCommendFromMachineLearnPai2%7Edefault-4.control
################################################

二、参数及输入文件说明
################################################
# 主要参数
# 一般情况下，我们只需要用到三个参数:
# -fasta 需要检索的蛋白序列的fasta文件
# -dir 存放Pfam-A数据库的目录
# -outfile 需要输出的文件名字
################################################
# 所有参数说明：
-dir                Pfam_data_file_dir   包含Pfam数据文件的目录[必须] 
-fasta              fasta_file   包含序列的输入文件名 [必须]
-e_seq              序列E-value阈值 [不指定则使用默认阈值] 
-e_dom              结构域E-value阈值 [不指定则使用默认阈值]
-b_seq              序列bit score阈值 [不指定则使用默认阈值]
-b_dom              结构域bit score阈值[不指定则使用默认阈值] 
-align              在结果中显示比对片段 [默认关闭] 
-as                 预测Pfam-A数据库匹配的active sites[默认关闭] 
-json [pretty]      输出结果使用JSON格式。例如指定值为[pretty]，则输出结果会使用"pretty" JSON格式输出 [默认关闭] 
-cpu                并行工作的CPU数目 [默认全部]
-translate [mode]   将输入序列视为DNA，并在搜索前使用6框翻译的方法进行转换。如果翻译模式[mode]被指定，则必须为"all"或者"orf"。"all"表示完整翻译，包括终止子并且不产生单独的ORFs；"orf"表示只翻译和报告长度大于20的ORFs。
如果使用了翻译参数而没有指定翻译模式，则默认使用"orf"模式。[默认关闭]
################################################

输入文件说明：
！注意：-fasta文件为如下所示的蛋白序列（只截取了前两列）
################################################
>S111001556  locus=S111-genome:1:1029:+
VHERSRLNPILTFDNLVTGKANQLARAAAVQVANNPGKSYNPLYLYGGVG
LGKTHLIHAIGNFMLMENPRARIRYIHAEQYVSDVVKAYQRKAFDDFKRY
YHSLDLLLIDDIQFFSGKNRTQEEFFYAFEALIANRAQVIITSDTYPKEI
TGIDDRLISRFDSGLTVAIEPPELEMRVAILMKKAQAENVTVPEEVAFFV
AKHLRSNVRELEGALRKILAYSNFHGKEITIEVTREALKDLLTVQNRQIS
VENIQKTCADFYNIKVADMYSKKRPANIARPRQIAMYLAKELTQKSLPEI
GELFGGRDHTTVLHAVRKIADERSKDAQLNHELHVLEQTLKG
>S111001557  locus=S111-genome:1307:2422:+
MQLVKTSRDNLLRPLQIVSGIVERRHTLPILANLLIRKNGERVSFLSTDI
EIQITTHADCGAGNGDIATTVAARKLVDILRAMPDGEVALTLNDKRMSVQ
SGKSRFALQTLAAEEFPTVAEANDFGARITLPQKTLKHLLAMVHFAMAQQ
DIRYYLNGMLLVVDGKQVMAVATDGHRLAYCGVETGEQPAGAGGRHEVII
PRKTILELQRLLEDVDDPVSVQLASNQVKFTFGNIELISKLVEGKFPDFQ
RVIPKGYRNSFTIDRAFLQQALQRTAILTTDKFKGVRCMLDTNVLKISST
NADQEEAQEELEIDYQGDALDIGFNVTYLLDVLANLKAEKVQVSLGDSNS
SALITLPDDDTFKYVVMPMRI
################################################

三、软件运行脚本
南农集群脚本
################################################
# 参考及脚本路径:/gss1/home/gfjiang/zhangyz/temp_Pfam/work.sh
source activate /gss1/home/gfjiang/anaconda3/envs/pfam
pfam_scan.pl -fasta /gss1/home/gfjiang/zhangyz/temp_Pfam/S111.faa  -dir /gss1/home/gfjiang/database/pfam/ -outfile /gss1/home/gfjiang/zhangyz/temp_Pfam/pfam_S111.faa -clan_overlap -as -cpu 16 -e_seq 1e-5 -e_dom 1e-5
################################################

pfam公司流程脚本(公司流程参考)
################################################
# 参考及脚本路径:/lustre/sdb/baopf/zyz/ykm/raw_protein/pfam.sh
cd /lustre/sdb/baopf/zyz/ykm/raw_protein
source activate /lustre/sdb/zengl/lib/anaconda3/envs/pfam
perl /lustre/sdb/zengl/bin/module/software/PfamScan/pfam_scan.pl -fasta RS-N.pep      -dir /lustre/sdb/zengl/bin/module/Anno/pfam/ -outfile pfam_RS-N.pep -clan_overlap -as -cpu 16 -e_seq 1e-5 -e_dom 1e-5
perl /lustre/sdb/zengl/bin/module/software/PfamScan/pfam_scan.pl -fasta S111.faa      -dir /lustre/sdb/zengl/bin/module/Anno/pfam/ -outfile pfam_S111.faa -clan_overlap -as -cpu 16 -e_seq 1e-5 -e_dom 1e-5
perl /lustre/sdb/zengl/bin/module/software/PfamScan/pfam_scan.pl -fasta YL-BAC-29.faa -dir /lustre/sdb/zengl/bin/module/Anno/pfam/ -outfile pfam_YL-BAC-29.faa -clan_overlap -as -cpu 16 -e_seq 1e-5 -e_dom 1e-5
perl /lustre/sdb/zengl/bin/module/software/PfamScan/pfam_scan.pl -fasta YL-ENT-31.faa -dir /lustre/sdb/zengl/bin/module/Anno/pfam/ -outfile pfam_YL-ENT-31.faa -clan_overlap -as -cpu 16 -e_seq 1e-5 -e_dom 1e-5
perl /lustre/sdb/zengl/bin/module/software/PfamScan/pfam_scan.pl -fasta YL-STE-01.faa -dir /lustre/sdb/zengl/bin/module/Anno/pfam/ -outfile pfam_YL-STE-01.faa -clan_overlap -as -cpu 16 -e_seq 1e-5 -e_dom 1e-5
################################################

三、补充说明
# Pfam_scan用与单菌的蛋白序列的蛋白结构域注释，而MetaCLADE用于分析宏基因组、宏转录组蛋白质结构域的注释流程，后续会添加该流程
# ①MetaCLADE基于Python，它直接接受短序列输入，用于定量分析宏基因组和宏转录组中的结构域（功能）；
# ②首先预测短序列可能编码的氨基酸序列，然后通过搜索与15000个Pfam结构域相关的超过两百万个概率模型来鉴定蛋白质结构域；
# ③其注释结果较InterProScan更为完整，与最先进的UProC和HMM-GRASPx方法相比，MetaCLADE也有所长；
# ④联合使用UProC和MetaCLADE，可以获得最佳的注释结果。