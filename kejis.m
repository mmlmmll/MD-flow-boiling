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
                        current_atom_data(:, :, i) = zeros(current_Natoms(i), 8); % 假设原子数据每行有8列信息（按实际调整）
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

% 根据原子类型设置原子质量（假设原子数据每行第2列代表原子类型，按实际调整）
atomMasses = zeros(size(atom_data, 1), 1);
atomTypeMassMap = containers.Map({2, 3, 4}, {12.0110, 18.9980, 1});
for t = 1:numTimeSteps
    % 获取当前时间步的原子数据（去除多余维度）
    currentAtomData = squeeze(atom_data(:, :, t)); 
    for atomIndex = 1:size(currentAtomData, 1)
        atomType = currentAtomData(atomIndex, 2);
        atomMasses(atomIndex) = atomTypeMassMap(atomType);
    end
    atom_data(:, :, t) = [currentAtomData(:, 1), atomMasses, currentAtomData(:, 3:end)];
end

% 用于存储不同时间步下每个分子的平动动能和转动动能，初始化时增加时间步维度
pdkeall = zeros(numTimeSteps, [], 1);
zdkeall = zeros(numTimeSteps, [], 1);
zke = zeros(numTimeSteps, [], 1);

% 计算每个时间步下各分子的平动动能和转动动能
for t = 1:numTimeSteps
    % 当前时间步的原子数据
    atomData = squeeze(atom_data(:, :, t)); 
    
    % 假设每个分子由5个原子组成
    numAtomsInMol = 5; 
    numMols = size(atomData, 1) / numAtomsInMol; % 分子数量，基于当前时间步原子数据计算
    if ~mod(size(atomData, 1), numAtomsInMol) == 0
        warning(['在时间步 ', num2str(t), ' 原子数量不能被每个分子的原子数整除，可能存在数据问题']);
    end
    
    % 用于存储当前时间步每个分子的平动动能和转动动能
    pdke = zeros(numMols, 1);
    zdke = zeros(numMols, 1);
    zmolke = zeros(numMols, 1);
    for molIndex = 1:numMols
        % 提取当前分子的原子数据
        startRow = (molIndex - 1) * numAtomsInMol + 1;
        endRow = molIndex * numAtomsInMol;
        molAtomsData = atomData(startRow:endRow, :);
        
        % 计算分子质心坐标
        totalMass = sum(atomMasses(startRow:endRow));
        xCenter = sum(atomMasses(startRow:endRow).* molAtomsData(:, 3)) / totalMass;
        yCenter = sum(atomMasses(startRow:endRow).* molAtomsData(:, 4)) / totalMass;
        zCenter = sum(atomMasses(startRow:endRow).* molAtomsData(:, 5)) / totalMass;
        
        % 计算分子质心速度
        vxCenter = sum(atomMasses(startRow:endRow).* molAtomsData(:, 6)) / totalMass;
        vyCenter = sum(atomMasses(startRow:endRow).* molAtomsData(:, 7)) / totalMass;
        vzCenter = sum(atomMasses(startRow:endRow).* molAtomsData(:, 8)) / totalMass;
        
        % 计算平动动能
        pdke(molIndex) = 0.5 * totalMass * (vxCenter^2 + vyCenter^2 + vzCenter^2);
          % 计算总动能
        zmolke(molIndex) = 0.5 * sum(atomMasses(startRow:endRow).* molAtomsData(:, 6).* molAtomsData(:, 6)+atomMasses(startRow:endRow).* molAtomsData(:, 7).* molAtomsData(:, 7)+atomMasses(startRow:endRow).* molAtomsData(:, 8).* molAtomsData(:, 8));
        % 计算转动惯量（考虑三维情况，更接近实际但仍简化示例）
        Ixx = 0;
        Iyy = 0;
        Izz = 0;
        omegaX1=0;
        omegaY1=0;
        omegaZ1=0;
        for atomIndex = 1:numAtomsInMol
            atomRow = startRow + atomIndex - 1;
            dx = molAtomsData(atomIndex, 3) - xCenter;
            dy = molAtomsData(atomIndex, 4) - yCenter;
            dz = molAtomsData(atomIndex, 5) - zCenter;
            Ixx = Ixx + atomMasses(atomRow) * (dy^2 + dz^2);
            Iyy = Iyy + atomMasses(atomRow) * (dx^2 + dz^2);
            Izz = Izz + atomMasses(atomRow) * (dx^2 + dy^2);
       
        
        % 这里简单假设角速度（示例中设为固定值，实际需准确计算）
        omegaX1 = omegaX1 + (molAtomsData(atomIndex, 6)-vxCenter)/sqrt(dx^2+dy^2+dz^2); 
        omegaY1 = omegaY1 + (molAtomsData(atomIndex, 7)-vxCenter)/sqrt(dx^2+dy^2+dz^2);
        omegaZ1 = omegaZ1 + (molAtomsData(atomIndex, 8)-vxCenter)/sqrt(dx^2+dy^2+dz^2);
        end
        omegaX=omegaX1/numAtomsInMol;
        omegaY=omegaY1/numAtomsInMol;
        omegaZ=omegaZ1/numAtomsInMol;
        omegaX2 =  (molAtomsData(:, 6)-vxCenter)/sqrt(dx^2+dy^2+dz^2); 
        omegaY2 =  (molAtomsData(:, 7)-vxCenter)/sqrt(dx^2+dy^2+dz^2);
        omegaZ2 =  (molAtomsData(:, 8)-vxCenter)/sqrt(dx^2+dy^2+dz^2);
        % 计算转动动能（对于非线性分子按对应公式修改这里的计算）
        zdke(molIndex) = 0.5 * (Ixx * omegaX^2 + Iyy * omegaY^2 + Izz * omegaZ^2);
    end
    
    % 将当前时间步每个分子的平动动能和转动动能存入总数组，注意按时间步存储
    pdkeall(t, 1:numMols, 1) = pdke;
    zdkeall(t, 1:numMols, 1) = zdke;
    zke(t, 1:numMols, 1)= zmolke;
end

% 计算每个时间步下平动动能和转动动能的平均值
avepdke = mean(pdkeall, 2);
avezdke = mean(zdkeall, 2);
sumpdke = sum(pdkeall, 2);
sumzdke = sum(zdkeall, 2);
sumzke = sum(zke,2);
% 输出结果到新的.data文件
outputFileName = 'results.data';
fid = fopen(outputFileName, 'w');
if fid == -1
    error('无法创建输出文件');
end

% 写入每个时间步下平动动能和转动动能的平均值数据
fprintf(fid, 'time pzke ppke sumzke sumpke sumzke\n');
for t = 1:numTimeSteps
    fprintf(fid, '%f %f %f %f %f %f\n',t, avezdke(t) ,avepdke(t),sumzdke(t) ,sumpdke(t) ,sumzke(t));
end

fclose(fid);

% 输出结果示例（可以根据需要进一步完善输出格式等）
disp('结果已输出到 results.data 文件');

% 可以进一步进行绘图等分析，比如绘制平动动能和转动动能随时间步的变化
figure;
plot(timestep, avepdke');
title('平动动能随时间步变化');
xlabel('时间步');
ylabel('平动动能');

figure;
plot(timestep, avezdke');
title('转动动能随时间步变化');
xlabel('时间步');
ylabel('转动动能');

figure;
plot(timestep, sumpdke');
title('平动动能随时间步变化');
xlabel('时间步');
ylabel('平动动能');

figure;
plot(timestep, sumzdke');
title('转动动能随时间步变化');
xlabel('时间步');
ylabel('转动动能');
figure;
plot(timestep, sumzke');
title('转动动能随时间步变化');
xlabel('时间步');
ylabel('转动动能');