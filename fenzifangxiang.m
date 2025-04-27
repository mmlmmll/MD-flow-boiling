clc;clear;
% 存储所有.dump文件名的单元格数组（需根据实际情况修改为包含所有文件名的单元格数组）
dumpFileNames = {'feiteng1.dump','feiteng2.dump','feiteng3.dump','feiteng4.dump','feiteng5.dump','feiteng6.dump','feiteng7.dump','feiteng8.dump','feiteng9.dump','feiteng10.dump'}; % 示例文件名，按实际替换

% 用于存储整合后的不同时间步下的相关数据（整合多个文件）
timestep = [];
Natoms = [];
x_bound = [];
y_bound = [];
z_bound = [];
atom_data = [];

% 循环读取每个.dump文件并整合数据
readFilesStartTime = cputime;
for fileIndex = 1:length(dumpFileNames)
    file = dumpFileNames{fileIndex};
    try
        dump = fopen(file, 'r');
    catch
        error(['文件 ' file ' 未找到！']);
    end
    
    % 临时存储当前文件不同时间步下的相关数据
    current_timestep = [];
    current_Natoms = [];
    current_x_bound = [];
    current_y_bound = [];
    current_z_bound = [];
    current_atom_data = [];
    
    i = 1;
    while feof(dump) == 0
        id = fgetl(dump);
        if (strncmpi(id, 'ITEM: TIMESTEP', numel('ITEM: TIMESTEP')))
            current_timestep(i) = str2num(fgetl(dump));
        else
            if (strncmpi(id, 'ITEM: NUMBER OF ATOMS', numel('ITEM: NUMBER OF ATOMS')))
                current_Natoms(i) = str2num(fgetl(dump));
            else
                if (strncmpi(id, 'ITEM: BOX BOUNDS', numel('ITEM: BOX BOUNDS')))
                    current_x_bound(i, :) = str2num(fgetl(dump));
                    current_y_bound(i, :) = str2num(fgetl(dump));
                    current_z_bound(i, :) = str2num(fgetl(dump));
                else
                    if (strcmpi(id(1:11), 'ITEM: ATOMS'))
                        current_atom_data(:, :, i) = zeros(current_Natoms(i), 8); 
                        for j = 1 : 1: current_Natoms(i)
                            current_atom_data(j, :, i) = str2num(fgetl(dump));
                        end
                        i = i + 1;
                    end
                end
            end
        end
    end
    fclose(dump);
    
    % 将当前文件的数据整合到总的数据中
    timestep =[timestep, current_timestep];
    Natoms = [Natoms, current_Natoms];
    x_bound = [x_bound, current_x_bound];
    y_bound = [y_bound, current_y_bound];
    z_bound = [z_bound, current_z_bound];
    atom_data = cat(3,atom_data, current_atom_data);
    readFilesEndTime = cputime;
    readFilesTime = (readFilesEndTime - readFilesStartTime) / 60;
    fprintf('读取所有文件所用时间是: %.1f mins.\n', readFilesTime);
end
% 去除原子类型为1的原子（假设原子数据每行第2列代表原子类型，按实际调整）
numTimeSteps = size(atom_data, 3);
newAtom_data = [];
for t = 1:numTimeSteps
    % 获取当前时间步的原子数据（去除多余维度）
    currentAtomData = squeeze(atom_data(:, :, t)); 
    % 筛选出原子类型不为1的原子数据
    filteredAtomData = currentAtomData(currentAtomData(:, 2) ~= 1, :); 
    
    % 将筛选后的数据添加到新的atom_data数组中
    newAtom_data = cat(3, newAtom_data, filteredAtomData);
end
atom_data = newAtom_data;

% 处理原子序号乱序问题，先对每个时间步下的原子按序号排序
for t = 1:numTimeSteps
    % 获取当前时间步的原子数据（去除多余维度）
    currentAtomData = squeeze(atom_data(:, :, t)); 
    % 按照原子序号（假设原子数据每行第1列是原子序号，按实际调整）进行排序
    sortedAtomData = sortrows(currentAtomData, 1); 
    atom_data(:, :, t) = sortedAtomData;
end
        % 计算夹角
  angle_degree=[ ];
  h=1;
for t = 1:10:numTimeSteps
    % 当前时间步的原子数据
    atomData = squeeze(atom_data(:, :, t)); 
    
    % 假设每个分子由5个原子组成
    numAtomsInMol = 5; 
    numMols = size(atomData, 1) / numAtomsInMol; % 分子数量，基于当前时间步原子数据计算
    if ~mod(size(atomData, 1), numAtomsInMol) == 0
        warning(['在时间步 ', num2str(t), ' 原子数量不能被每个分子的原子数整除，可能存在数据问题']);
    end
    
     for molIndex = 1:numMols
        % 提取当前分子的原子数据
        startRow = (molIndex - 1) * numAtomsInMol + 1;
        endRow = molIndex * numAtomsInMol;
        molAtomsData = atomData(startRow:endRow, :);
     
            A = [molAtomsData(1,3), molAtomsData(1,4), molAtomsData(1,5)];
            B = [molAtomsData(2,3), molAtomsData(2,4), molAtomsData(2,5)];
            C = [molAtomsData(3,3), molAtomsData(3,4), molAtomsData(3,5)];
            % 计算向量 AB 和 AC
            AB = B - A;
            AC = C - A;
        % 计算面的法向量 n，通过叉乘得到
            ab=AB/norm(AB);ac=AC/norm(AC);
            N=ab+ac;
            % 对法向量进行归一化处理
            n = N / norm(N);   
            

            % z 轴的单位向量 k
            k = [0, 0, 1];
            % 计算两个向量的点积
            dot_product = dot(n, k);
            % 计算夹角，使用 acos 函数，结果为弧度
            angle = acos(dot_product);
            % 将弧度转换为角度
           angeles = rad2deg(angle);
            angle_degree(molIndex,h)=angeles;            
     end
      h=h+1;
end
 % 存储结果到文件
dlmwrite('matrix_data.txt', angle_degree, 'delimiter', ' ');