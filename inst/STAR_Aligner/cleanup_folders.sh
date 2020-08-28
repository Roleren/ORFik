usage(){
cat << EOF
usage: $0 options

script to clean up after STAR run

OPTIONS:
	input: path to folder of STAR output files
EOF
}

out_dir=$1
echo "cleaning up $out_dir"
if [ ! -d ${out_dir} ]; then
        echo "$out_dir does not exist!"
        exit 1
fi

if [ ! -d "${out_dir}/aligned" ]; then
        echo "${out_dir}/aligned does not exist!"
        exit 1
fi

# clean up

# temp files

# Logs to seperate folder
function cleanup_all()
{

	if [ -d $1 ]; then
		if [ $(ls ${1} | wc -l)  -gt 0 ]; then
			find $1 -maxdepth 1 | grep '_STARtmp' | xargs rm -r
			if [ ! -d $1/LOGS ]; then
				mkdir $1/LOGS
			fi
			find $1 -maxdepth 1 | grep '_Log\.' | xargs -I {} mv {} $1/LOGS/
			find $1 -maxdepth 1 | grep '\.out.tab' | xargs -I {} mv {} $1/LOGS/
		fi
	fi
}

cleanup_all ${out_dir}/aligned
cleanup_all ${out_dir}/ncRNA_depletion
cleanup_all ${out_dir}/phix_depletion
cleanup_all ${out_dir}/rRNA_depletion
cleanup_all ${out_dir}/tRNA_depletion



