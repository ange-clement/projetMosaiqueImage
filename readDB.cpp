#include <iostream>
#include <sqlite3.h>
#include <string>
  
using namespace std;
  
static int callback(void* data, int argc, char** argv, char** azColName)
{
    int i;
    fprintf(stderr, "%s: ", (const char*)data);
  
    for (i = 0; i < argc; i++) {
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
  
    printf("\n");
    return 0;
}
  
int main(int argc, char** argv)
{
    sqlite3* DB;
    char* messaggeError;
    int exit = sqlite3_open("example.db", &DB);
    string query = "SELECT * FROM IMAGE;";
  
    cout << "STATE OF TABLE BEFORE INSERT" << endl;
  
    sqlite3_exec(DB, query.c_str(), callback, NULL, NULL);
    return(0);
  }