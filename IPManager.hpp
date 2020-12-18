#ifndef IPMANAGER_H__
#define IPMANAGER_H__

#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
using namespace std;
using namespace boost;


class IPManager 
{
public:
	int machines;
	string A_Ip[100];
	string B_Ip[100];
	string C_Ip[100];
	string D_Ip[100];
	string P_Ip[4][100];

	
	unordered_map <string, string> machineIps;

	IPManager(int machines){
		loadIPs(machines);
	}

	~IPManager(){
	}

	void loadIPs(int machines) {
		ifstream file ("./machine_spec/" + to_string(machines) + ".txt"); 
		if( !file ) {
			cout << "Error opening input file" << endl ;
		}
   		string line = "";
   		while (file.peek() != EOF) {
   			getline (file , line);
   			while (starts_with(line, "#")) {
				getline(file, line);
			}
			if (starts_with(line, "processor_spec")) {
				break;
			}
			// cout << line << "\n";
			vector<string> token;
			split(token, line, is_any_of(" "), token_compress_on);
			machineIps.emplace(token[0], token[1]);
     	}
     	getline(file, line);
     	while (starts_with(line, "#")) {
			getline(file, line);
		}
		// cout << line << "\n";
     	setIps(machines, line);
		getline(file, line);
     	while (starts_with(line, "#")) {
			getline(file, line);
		}
		// cout << line << "\n";
     	setIps(machines, line);
     	getline(file, line);
     	while (starts_with(line, "#")) {
			getline(file, line);
		}
		// cout << line << "\n";
     	setIps(machines, line);
     	getline(file, line);
     	while (starts_with(line, "#")) {
			getline(file, line);
		}
		// cout << line << "\n";
     	setIps(machines, line);
		file.close();
	}

	// void checkComment(ifstream file, string next) {
	// 	getline(file, next);
	// 	// while (next.starts_with("#")) {
	// 	// 	getline(file, next);
	// 	// }
	// 	return next;
	// }

	void setIps(int machines, string next) {
		if (starts_with(next, "AliceMachines")) {
			vector<string> tok;
			split(tok, next, is_any_of(" "), token_compress_on);
			if (tok.size() != (machines+1)) {
				cout << "You have specified " << sizeof(tok)-1 << "  machines for Alice instead of " << machines << "\n";
			}
			for (int i = 1; i < tok.size(); i++) {
				P_Ip[0][i-1] = machineIps[tok[i]];
				// A_Ip[i-1] = machineIps[tok[i]];
				// cout << A_Ip[i-1] << " " <<P_Ip[0][i-1] << "\n";
				
			}
		} 
		else if (starts_with(next, "BobMachines")) {
			vector<string> tok;
			split(tok, next, is_any_of(" "), token_compress_on);
			if (tok.size() != (machines+1)) {
				cout << "You have specified " << sizeof(tok)-1 << "  machines for Bob instead of " << machines << "\n";
			}
			for (int i = 1; i < tok.size(); i++) {
				P_Ip[1][i-1] = machineIps[tok[i]];
				// B_Ip[i-1] = machineIps[tok[i]];
				// cout << B_Ip[i-1] << " " <<P_Ip[1][i-1] << "\n";
			}
		}
		else if (starts_with(next, "CharlieMachines")) {
			vector<string> tok;
			split(tok, next, is_any_of(" "), token_compress_on);
			if (tok.size() != (machines+1)) {
				cout << "You have specified " << sizeof(tok)-1 << "  machines for Charlie instead of " << machines << "\n";
			}
			for (int i = 1; i < tok.size(); i++) {
				P_Ip[2][i-1] = machineIps[tok[i]];
				// C_Ip[i-1] = machineIps[tok[i]];
				// cout << C_Ip[i-1] << " " <<P_Ip[2][i-1] << "\n";
			}
		}
// 		else if (starts_with(next, "DavidMachines")) {
// 			vector<string> tok;
// 			split(tok, next, is_any_of(" "), token_compress_on);
// 			if (tok.size() != (machines+1)) {
// 				cout << "You have specified " << sizeof(tok)-1 << "  machines for David instead of " << machines << "\n";
// 			}
// 			for (int i = 1; i < tok.size(); i++) {
// 				P_Ip[3][i-1] = machineIps[tok[i]];
// 				// D_Ip[i-1] = machineIps[tok[i]];
// 				// cout << D_Ip[i-1] << " " <<P_Ip[3][i-1] << "\n";
// 			}
// 		}

	}
};


#endif