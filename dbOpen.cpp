#include <iostream>
#include <sqlite3.h>
#include <string>
#include <filesystem>
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    if(argc<2){
        std::cout<<"Argument : databaseName"<<std::endl;
        exit(0);
    }

    std::string dbName(argv[1]);
    dbName+=".db";

    sqlite3* DB;    
    char* messaggeError;
    std::string sql = "CREATE TABLE IMAGE("
                      "ID INT PRIMARY KEY     NOT NULL, "
                      "NAME           TEXT    NOT NULL, "
                      "MEAN           REAL    NOT NULL);";

    int exit = sqlite3_open(dbName.c_str(), &DB);
    exit = sqlite3_exec(DB, sql.c_str(), NULL, 0, &messaggeError);
    if (exit) {
        std::cerr << "Error open DB " << sqlite3_errmsg(DB) << std::endl;
        return (-1);
    }
    else{
        std::cout << "Opened Database Successfully!" << std::endl;

        std::string path = "BOWS/";
        int i=0;
        for (const auto & entry : fs::directory_iterator(path)){
            std::string filename = entry.path();

            std::string sql("INSERT INTO IMAGE VALUES(");
            sql = sql+std::to_string(i)+",'"+filename+"', 0 );";
            exit = sqlite3_exec(DB, sql.c_str(), NULL, 0, &messaggeError);
            std::cout<<sql<<std::endl;
            if (exit != SQLITE_OK) {
                std::cerr << "Error Insert" <<messaggeError<< std::endl;
            }
            i++;
        }
    }
    sqlite3_close(DB);
    return (0);
}