clc;clear;
for n=1:10
     original_file_name = sprintf('feiteng%d.dump', n);
file ="original_file_name";

try
    dump = fopen(original_file_name,'r');
catch
    error('Dumpfile not found!');
end

i=1;
while feof(dump) == 0
    id = fgetl(dump);
     if (strncmpi(id,'ITEM: TIMESTEP',numel('ITEM: TIMESTEP')))
            timestep(i) = str2num(fgetl(dump));
     else
     if (strncmpi(id,'ITEM: NUMBER OF ATOMS',numel('ITEM: NUMBER OF ATOMS')))
            Natoms(i) = str2num(fgetl(dump));
     else
      if (strncmpi(id,'ITEM: BOX BOUNDS',numel('ITEM: BOX BOUNDS')))
            x_bound(i,:) = str2num(fgetl(dump));
            y_bound(i,:) = str2num(fgetl(dump));
            z_bound(i,:) = str2num(fgetl(dump));
      else
       if (strcmpi(id(1:11),'ITEM: ATOMS'))
            for j = 1 : 1: Natoms
                atom_data(j,:,i) = str2num(fgetl(dump));
            end
            i=i+1;
       end
      end 
     end
   end
end


%======================================%
cutoff   = 5.3;
 processed_file_name = sprintf('%d.dump', n);
[fid,message]  = fopen(processed_file_name,"wt");
%======================================%
all_frame = size(atom_data,3);
xl = x_bound(1,2)-x_bound(1,1);
yl = y_bound(1,2)-y_bound(1,1);
zl = z_bound(1,2)-z_bound(1,1);

Ntotal   = length(atom_data);
nerghbor_num = zeros(length(atom_data),size(atom_data,3));
gas_num = zeros(length(atom_data),size(atom_data,3));
distance  = zeros(length(atom_data),1);


for frame = 1:all_frame
    
    now_frame = atom_data(:,:,frame);
    
    ID        = now_frame(:,1);
    TYPE      = now_frame(:,2);
    XYZ       = now_frame(:,3:5);
    cluster_tag = zeros(Ntotal,1);
    type = TYPE;

    %=====================================================%
    parfor i = 1:Ntotal   

        for j = 1:Ntotal
            if(i==j||TYPE(i)~=2||type(j)~=2)
                continue;
            end
            
            dx = (XYZ(j,1)-XYZ(i,1));
            dy = (XYZ(j,2)-XYZ(i,2));            
            dz = (XYZ(j,3)-XYZ(i,3));
           if(dx>xl/2)
        dx = dx-xl;
            elseif (dx<-xl/2)
        dx = dx+xl;
           end
    

           if(dy>yl/2)
        dy = dy-yl;
    elseif (dy<-yl/2)
        dy = dy+yl;
           end
        
            
            distance(i) = sqrt(dx*dx+dy*dy+dz*dz);
            
            if(distance(i)<cutoff&&TYPE(i)==2)
                nerghbor_num(i,frame)=nerghbor_num(i,frame)+1;    
            end
            
        end
        
        if(nerghbor_num(i,frame)<4&&TYPE(i)==2)
            TYPE(i)=5;   
        end
        
    end % find nerghbor 
    %=====================================================%
    
     
 %写出dump轨迹    
 write_data      = zeros(Ntotal,6); 

 write_data(:,1)   = ID;       
 write_data(:,2)   = TYPE;
 write_data(:,3:5) = XYZ;  
% %------------------------------------------------%  
    title_first    = "ITEM: TIMESTEP";
    title_number   = "ITEM: NUMBER OF ATOMS";
    title_boundary = "ITEM: BOX BOUNDS pp pp pp";
    title_all      = "ITEM: ATOMS id type x y z";
         
    fprintf(fid,'%s\n',title_first);
    fprintf(fid,'%d\n',timestep(frame));
    
    fprintf(fid,'%s\n',title_number);
    fprintf(fid,'%d\n',Natoms(1));
    
    fprintf(fid,'%s\n',title_boundary);
    fprintf(fid,'%d %d\n',[x_bound(1,1),x_bound(1,2)]);
    fprintf(fid,'%d %d\n',[y_bound(1,1),y_bound(1,2)]);
    fprintf(fid,'%d %d\n',[z_bound(1,1),z_bound(1,2)]);     
        
    fprintf(fid,'%s\n',title_all);
    [r,c]=size(write_data(:,:));         
        for i=1:r
            for j=1:c
                fprintf(fid,'%f\t',write_data(i,j));
            end
                fprintf(fid,'\r\n');
        end 
%------------------------------------------------%      
fprintf('Now the frame is: %.1f.\n',frame);
disp("-------------------");
end % loop all frame 

disp("-------------------");
disp("----ALL DONE!!!%d----");
disp("-------------------");
end

