clc;clear;
% �洢����.dump�ļ����ĵ�Ԫ�����飨�����ʵ������޸�Ϊ���������ļ����ĵ�Ԫ�����飩
dumpFileNames = {'feiteng1.dump','feiteng2.dump','feiteng3.dump','feiteng4.dump','feiteng5.dump','feiteng6.dump','feiteng7.dump','feiteng8.dump','feiteng9.dump','feiteng10.dump'}; % ʾ���ļ�������ʵ���滻

% ���ڴ洢���Ϻ�Ĳ�ͬʱ�䲽�µ�������ݣ����϶���ļ���
timestep = [];
Natoms = [];
x_bound = [];
y_bound = [];
z_bound = [];
atom_data = [];

% ѭ����ȡÿ��.dump�ļ�����������
readFilesStartTime = cputime;
for fileIndex = 1:length(dumpFileNames)
    file = dumpFileNames{fileIndex};
    try
        dump = fopen(file, 'r');
    catch
        error(['�ļ� ' file ' δ�ҵ���']);
    end
    
    % ��ʱ�洢��ǰ�ļ���ͬʱ�䲽�µ��������
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
    
    % ����ǰ�ļ����������ϵ��ܵ�������
    timestep =[timestep, current_timestep];
    Natoms = [Natoms, current_Natoms];
    x_bound = [x_bound, current_x_bound];
    y_bound = [y_bound, current_y_bound];
    z_bound = [z_bound, current_z_bound];
    atom_data = cat(3,atom_data, current_atom_data);
    readFilesEndTime = cputime;
    readFilesTime = (readFilesEndTime - readFilesStartTime) / 60;
    fprintf('��ȡ�����ļ�����ʱ����: %.1f mins.\n', readFilesTime);
end
% ȥ��ԭ������Ϊ1��ԭ�ӣ�����ԭ������ÿ�е�2�д���ԭ�����ͣ���ʵ�ʵ�����
numTimeSteps = size(atom_data, 3);
newAtom_data = [];
for t = 1:numTimeSteps
    % ��ȡ��ǰʱ�䲽��ԭ�����ݣ�ȥ������ά�ȣ�
    currentAtomData = squeeze(atom_data(:, :, t)); 
    % ɸѡ��ԭ�����Ͳ�Ϊ1��ԭ������
    filteredAtomData = currentAtomData(currentAtomData(:, 2) ~= 1, :); 
    
    % ��ɸѡ���������ӵ��µ�atom_data������
    newAtom_data = cat(3, newAtom_data, filteredAtomData);
end
atom_data = newAtom_data;

% ����ԭ������������⣬�ȶ�ÿ��ʱ�䲽�µ�ԭ�Ӱ��������
for t = 1:numTimeSteps
    % ��ȡ��ǰʱ�䲽��ԭ�����ݣ�ȥ������ά�ȣ�
    currentAtomData = squeeze(atom_data(:, :, t)); 
    % ����ԭ����ţ�����ԭ������ÿ�е�1����ԭ����ţ���ʵ�ʵ�������������
    sortedAtomData = sortrows(currentAtomData, 1); 
    atom_data(:, :, t) = sortedAtomData;
end
        % ����н�
  angle_degree=[ ];
  h=1;
for t = 1:10:numTimeSteps
    % ��ǰʱ�䲽��ԭ������
    atomData = squeeze(atom_data(:, :, t)); 
    
    % ����ÿ��������5��ԭ�����
    numAtomsInMol = 5; 
    numMols = size(atomData, 1) / numAtomsInMol; % �������������ڵ�ǰʱ�䲽ԭ�����ݼ���
    if ~mod(size(atomData, 1), numAtomsInMol) == 0
        warning(['��ʱ�䲽 ', num2str(t), ' ԭ���������ܱ�ÿ�����ӵ�ԭ�������������ܴ�����������']);
    end
    
     for molIndex = 1:numMols
        % ��ȡ��ǰ���ӵ�ԭ������
        startRow = (molIndex - 1) * numAtomsInMol + 1;
        endRow = molIndex * numAtomsInMol;
        molAtomsData = atomData(startRow:endRow, :);
     
            A = [molAtomsData(1,3), molAtomsData(1,4), molAtomsData(1,5)];
            B = [molAtomsData(2,3), molAtomsData(2,4), molAtomsData(2,5)];
            C = [molAtomsData(3,3), molAtomsData(3,4), molAtomsData(3,5)];
            % �������� AB �� AC
            AB = B - A;
            AC = C - A;
        % ������ķ����� n��ͨ����˵õ�
            ab=AB/norm(AB);ac=AC/norm(AC);
            N=ab+ac;
            % �Է��������й�һ������
            n = N / norm(N);   
            

            % z ��ĵ�λ���� k
            k = [0, 0, 1];
            % �������������ĵ��
            dot_product = dot(n, k);
            % ����нǣ�ʹ�� acos ���������Ϊ����
            angle = acos(dot_product);
            % ������ת��Ϊ�Ƕ�
           angeles = rad2deg(angle);
            angle_degree(molIndex,h)=angeles;            
     end
      h=h+1;
end
 % �洢������ļ�
dlmwrite('matrix_data.txt', angle_degree, 'delimiter', ' ');