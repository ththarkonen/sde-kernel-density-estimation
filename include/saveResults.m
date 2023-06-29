function [] = saveResults( modelString, theta, y)

    timeStamp = datetime( "now", "Format", "yyyy_MM_dd_hh_mm_ss" );
    saveDirPath = "./results/" + modelString + "-" + string( timeStamp );

    thetaFilePath = saveDirPath + "/" + modelString + "-theta-results.mat";
    modelFilePath = saveDirPath + "/" + modelString + "-data.mat";

    mkdir( saveDirPath );
    
    save( thetaFilePath, "theta")
    save( modelFilePath, "y")
end