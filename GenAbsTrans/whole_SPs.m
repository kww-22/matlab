
humSPif = [];
humSPof = [];
faSPif = [];
faSPof = [];

mhumSPif = [];
mhumSPof = [];
mfaSPif = [];
mfaSPof = [];

for i = 1:numFiles
    data = extractData(fileNames.fileNames{i},'text',9);
    
    events = find(data.VEM_0 == 1);
    data = data(events(1):events(end),:);
    
    humSPif = [humSPif; max(data.HumSp)];
    faSPif = [faSPif; max(data.FaSp)];
    humSPof = [humSPof; min(data.HumSp)];
    faSPof = [faSPof; min(data.FaSp)];

    mhumSPif = [mhumSPif; max(data.HumSp)/max(data.BodyMass)];
    mfaSPif = [mfaSPif; max(data.FaSp)/max(data.BodyMass)];
    mhumSPof = [mhumSPof; min(data.HumSp)/max(data.BodyMass)];
    mfaSPof = [mfaSPof; min(data.FaSp)/max(data.BodyMass)];

end