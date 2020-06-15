#usage
#gawk -f count.awk hg19.fa
BEGIN{
    window = 10;
}

/>/{
    title = $0;
    line = 0;
}

!/>/{
    line ++;
    data = data $0;
    for(i=1;i<=length(data)-window;i++){
        hash = substr(data,i,window);
        #print line, hash;
        count[hash] ++;
    }
    data = substr(data,length(data)-window+1); #next char
    #print line,data;
}

END{
    count[data]++;
    for(tag in count){
        print tag,count[hash]
    }
}
