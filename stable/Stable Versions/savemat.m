function ty = savemat( fnaam, A )
%savemat(name,matrix)
%last change 30-12-99

[a,b]=size(A);
if b>1
    form='%e\t ';
    for i=1:b-2
        form=[form,'%e\t '];
    end
    form=[form,'%e\r\n'];
else
    form=['%e\r\n'];
end
fid=fopen(fnaam,'w');
if fid~=-1
    fprintf(fid,form,A');
    fclose(fid);
end
ty=fid;   