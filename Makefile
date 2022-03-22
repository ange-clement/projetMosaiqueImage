reader:
	c++ readDB.cpp -l sqlite3 -o read -std=c++17
creator:
	c++ dbOpen.cpp -l sqlite3 -o create -std=c++17

