/*
 * Copyright (C) 2010-2012 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CONFIG_FILE_HH
#define _CONFIG_FILE_HH

#include <map>
#include <vector>

#include "WarningConsolidation.hh"

using namespace std;

typedef map<string, double> member_nodes_t;

class Hook
{
public:
  Hook(string &global_init_file_arg);
  ~Hook();
private:
  map<string, string> hooks;
public:
  inline map<string, string>get_hooks() { return hooks; };
};

class SlaveNode
{
  friend class ConfigFile;
public:
  SlaveNode(string &computerName_arg, string port_arg, int minCpuNbr_arg, int maxCpuNbr_arg, string &userName_arg,
            string &password_arg, string &remoteDrive_arg, string &remoteDirectory_arg,
            string &dynarePath_arg, string &matlabOctavePath_arg, bool singleCompThread_arg,
            string &operatingSystem_arg);
  ~SlaveNode();

protected:
  const string computerName;
  const string port;
  int minCpuNbr;
  int maxCpuNbr;
  const string userName;
  const string password;
  const string remoteDrive;
  const string remoteDirectory;
  const string dynarePath;
  const string matlabOctavePath;
  const bool singleCompThread;
  const string operatingSystem;
};

class Cluster
{
  friend class ConfigFile;
public:
  Cluster(member_nodes_t member_nodes_arg);
  ~Cluster();

protected:
  member_nodes_t member_nodes;
};

//! The abstract representation of a "config" file
class ConfigFile
{
public:
  ConfigFile(bool parallel_arg, bool parallel_test_arg, bool parallel_slave_open_mode_arg, const string &cluster_name, int dct_arg);
  ~ConfigFile();

private:
  const bool parallel;
  const bool parallel_test;
  const bool parallel_slave_open_mode;
  const int dct;	
  const string cluster_name;
  string firstClusterName;
  //! Hooks
  vector<Hook *> hooks;
  //! Cluster Table
  map<string, Cluster *> clusters;
  //! Node Map
  map<string, SlaveNode *> slave_nodes;
  //! Add Hooks
  void addHooksConfFileElement(string &global_init_file);
  //! Add a SlaveNode or a Cluster object
  void addParallelConfFileElement(bool inNode, bool inCluster, member_nodes_t member_nodes, string &name,
                                  string &computerName, string port, int minCpuNbr, int maxCpuNbr, string &userName,
                                  string &password, string &remoteDrive, string &remoteDirectory,
                                  string &dynarePath, string &matlabOctavePath, bool singleCompThread,
                                  string &operatingSystem);
public:
  //! Parse config file
  void getConfigFileInfo(const string &parallel_config_file);
  //! Check Pass
  void checkPass(WarningConsolidation &warnings) const;
  //! Check Pass
  void transformPass();
  //! Write any hooks
  void writeHooks(ostream &output) const;
  //! Create options_.parallel structure, write options
  void writeCluster(ostream &output) const;
  //! Close slave nodes if needed
  void writeEndParallel(ostream &output) const;
};

#endif // ! CONFIG_FILE_HH
