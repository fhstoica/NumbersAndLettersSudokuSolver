//Compile under cygwin: g++ -Wall -std=gnu++0x -o NumberOrLetterSudokuSolver.srl NumberOrLetterSudokuSolver.cpp -I/cygdrive/c/cygwin/usr/include/boost

#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include <math.h>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

template <typename T> class Node{
public:
  typedef typename boost::shared_ptr< std::set<T> > SetPtrType;
private:
  unsigned int row;
  unsigned int col;
  bool empty;
  T value;
  static std::set<T> allValues; 
  SetPtrType possibleValues;
  SetPtrType rowSet;
  SetPtrType colSet;
  SetPtrType squareSet;
  void genZeroVal(T& val);
  void genPossibleValues(SetPtrType& S);
  T zero_value;
public:
  Node(unsigned int& row_pos, unsigned int& col_pos) : row(row_pos), col(col_pos) {
    rowSet = SetPtrType(new typename std::set<T>());
    colSet = SetPtrType(new typename std::set<T>());
    squareSet = SetPtrType(new typename std::set<T>());
    possibleValues = SetPtrType(new typename std::set<T>());
    empty = true;
    genZeroVal(zero_value); //Empty squares are marked with 0 for int and "_" for string.
    value = zero_value;
  }
  ~Node(){
  }
  SetPtrType& getRowSet(void);
  SetPtrType& getColSet(void);
  SetPtrType& getSquareSet(void);
  SetPtrType& getPossibleValues(void);
  void updatePossibleValues();
  void setValue(T& val);
  T getValue(void);
  bool isEmpty(void){return(empty);};
};

/************* Template specialization **************/
template<> void Node<unsigned int>::genZeroVal(unsigned int& val){
  val = static_cast<unsigned int>(0);
}
template<> void Node<std::string>::genZeroVal(std::string& val){
  val = std::string("_");
}
/****************************************************/
template<typename T> typename Node<T>::SetPtrType& Node<T>::getRowSet(void){
  return(rowSet);
}

template<typename T> typename Node<T>::SetPtrType& Node<T>::getColSet(void){
  return(colSet);
}
template<typename T> typename Node<T>::SetPtrType& Node<T>::getSquareSet(void){
  return(squareSet);
}
template<typename T> typename Node<T>::SetPtrType& Node<T>::getPossibleValues(void){
  return(possibleValues);
}
/***** Template specialization for static members; Note that this requires -std=gnu++0x with GCC 4.5 series *****/
template<> std::set<unsigned int> Node<unsigned int>::allValues({1, 2, 3, 4, 5, 6, 7 ,8, 9});
template<> std::set<std::string>  Node<std::string>::allValues({std::string("A"), std::string("B"), std::string("C"), std::string("D"),
      std::string("E"), std::string("F"), std::string("G"), std::string("H"), 
      std::string("I")});
/***********************************************************************************************************/

template<typename T> void Node<T>::updatePossibleValues(void){
  for(typename std::set<T>::iterator it = allValues.begin(); it != allValues.end(); ++it){
    possibleValues->insert(*it);
  }
  /* std::set_union() doesn't work with std::set ??? Apparently only with vectors. Strange
  SetPtrType rowColSet = SetPtrType(new typename std::set<T>());
  std::set_union(rowSet->begin(), rowSet->end(), colSet->begin(), colSet->end(), rowColSet->begin());

  SetPtrType rowColSqSet = SetPtrType(new typename std::set<T>());
  std::set_union(rowColSet->begin(), rowColSet->end(), squareSet->begin(), squareSet->end(), rowColSqSet->begin());
  */
  for(typename SetPtrType::element_type::iterator it = rowSet->begin(); it != rowSet->end(); ++it){
    if(possibleValues->find(*it) != possibleValues->end()){
      possibleValues->erase(*it);
    }
  }
  for(typename SetPtrType::element_type::iterator it = colSet->begin(); it != colSet->end(); ++it){
    if(possibleValues->find(*it) != possibleValues->end()){
      possibleValues->erase(*it);
    }
  }
  for(typename SetPtrType::element_type::iterator it = squareSet->begin(); it != squareSet->end(); ++it){
    if(possibleValues->find(*it) != possibleValues->end()){
      possibleValues->erase(*it);
    }
  }
}
template<typename T> void Node<T>::setValue(T& val){
  value = val;
  if(zero_value == value){
    empty = true;
  }
  else{
    empty = false;
  }
}
template<typename T> T Node<T>::getValue(void){
  return(value);
}

template<typename T> class Table{
private:
  const unsigned int N;
  const unsigned int max_depth;
  unsigned int sqSize;
  static bool solution_found;
  static bool valid_table;
  T zero_value;
  boost::shared_ptr<std::vector< std::vector<Node<T> > > > NodeTable;
  std::vector<typename Node<T>::SetPtrType> rowSets;
  std::vector<typename Node<T>::SetPtrType> colSets;
  std::vector< std::vector<typename Node<T>::SetPtrType > > squareSets;
  std::vector<unsigned int> rowSizeVec;
  std::vector<unsigned int> colSizeVec;
  std::vector< std::vector< unsigned int > > squareSizeVec;
  void genZeroVal(T& val);
  bool checkIfSolved(std::vector< std::vector<Node<T> > >& recNodeTable);
  void recursiveGenValueComb(std::vector<std::pair<unsigned int, unsigned int> > nodeVector,
			     std::vector<std::pair<unsigned int, unsigned int> > visitedNodeVector,
			     std::vector<T>& valueVector,
			     std::vector<std::vector<T> >& valueCombinations);
  void recursiveSolve(unsigned int depth);
  bool solveSinglePossValues( std::vector<std::pair<unsigned int, unsigned int> >& singleValuesSolved );
  void reverseSinglePossValues( std::vector<std::pair<unsigned int, unsigned int> >& singleValuesSolved );
public:
  Table(unsigned int size);
  ~Table();
  void updateNodeTable(std::vector<std::string> string_table);
  void updateSets();
  void solvePuzzle(void);
  std::string tableToString(void);
};
/**** Initialize static data *****/
template<typename T> bool Table<T>::solution_found = false;
template<typename T> bool Table<T>::valid_table    = true;
/*********************************/
template<typename T> Table<T>::Table(unsigned int size) : N(size), max_depth(N*N) {
  NodeTable = boost::shared_ptr<std::vector< std::vector<Node<T> > > >(new std::vector< std::vector<Node<T> > >());
  genZeroVal(zero_value);
  for(unsigned int i = 0; i < N; ++i){
    rowSizeVec.push_back(0);
    colSizeVec.push_back(0);
    rowSets.push_back(typename Node<T>::SetPtrType(new typename std::set<T>()));
    colSets.push_back(typename Node<T>::SetPtrType(new typename std::set<T>()));
  }
  sqSize = static_cast<unsigned int>(sqrt(N));
  if(static_cast<unsigned int>(N/sqSize) != sqSize){
    /*N must be a perfect square*/
    std::cerr << " N is not a perfect square : N = " << N << " sqSize = " << sqSize << std::endl;
    exit(1);
  }
  else{
    for(unsigned int j = 0; j < sqSize; ++j){
      std::vector<typename Node<T>::SetPtrType> tempSets;
      std::vector<unsigned int> tempVec;
      for(unsigned int k = 0; k < sqSize; ++k){
	tempSets.push_back(typename Node<T>::SetPtrType(new typename std::set<T>()));
	tempVec.push_back(0);
      }
      squareSets.push_back(tempSets);
      squareSizeVec.push_back(tempVec);
    }
  }
}

template<typename T> Table<T>::~Table(){
  rowSets.clear();
  colSets.clear();
  squareSets.clear();
}

/************* Template specialization **************/
template<> void Table<unsigned int>::genZeroVal(unsigned int& val){
  val = static_cast<unsigned int>(0);
}
template<> void Table<std::string>::genZeroVal(std::string& val){
  val = std::string("_");
}

template<> void Table<unsigned int>::updateNodeTable(std::vector<std::string> string_table){
  if(N != string_table.size()){
    std::cerr << "ERROR: N = " << N << " string_table.size() = " << string_table.size() << " exiting " << std::endl;
    exit(-1);
  }
  std::vector< Node<unsigned int> > nodeRow;
  for(unsigned int row = 0; row < N; ++row){
    if(N != string_table[row].size()){
      std::cerr << "ERROR: N = " << N 
		<< " string_table[" << row << "].size() = " 
		<< string_table.size() << " exiting " << std::endl;
      exit(-1);
    }
    nodeRow.clear();
    for(unsigned int col = 0; col < N; ++col){
      Node<unsigned int> tmp_node(row, col);
      char cc = string_table[row].substr(col, 1).c_str()[0] - (char)('0');
      unsigned int tmp_val = static_cast<unsigned int>(cc);
      tmp_node.setValue(tmp_val);
      nodeRow.push_back(tmp_node);
    }
    NodeTable->push_back(nodeRow);
  }
}

template<> void Table<std::string>::updateNodeTable(std::vector<std::string> string_table){
  if(N != string_table.size()){
    std::cerr << "ERROR: N = " << N << " string_table.size() = " << string_table.size() << " exiting " << std::endl;
    exit(-1);
  }
  std::vector< Node<std::string> > nodeRow;
  for(unsigned int row = 0; row < N; ++row){
    if(N != string_table[row].size()){
      std::cerr << "ERROR: N = " << N 
		<< " string_table[" << row << "].size() = " 
		<< string_table.size() << " exiting " << std::endl;
      exit(-1);
    }
    nodeRow.clear();
    for(unsigned int col = 0; col < N; ++col){
      Node<std::string> tmp_node(row, col);
      std::string cc = string_table[row].substr(col, 1);
      std::string tmp_val(cc);
      tmp_node.setValue(tmp_val);
      nodeRow.push_back(tmp_node);
    }
    NodeTable->push_back(nodeRow);
  }
}
/****************************************************/

template<typename T> void Table<T>::updateSets(){

  /* Also check if the configuration is a valid table, i.e. no repeated values. */
  rowSets.clear();
  colSets.clear();
  squareSets.clear();

  rowSizeVec.clear();
  colSizeVec.clear();
  squareSizeVec.clear();

  valid_table = true;

  for(unsigned int i = 0; i < N; ++i){
    rowSets.push_back(typename Node<T>::SetPtrType(new typename std::set<T>()));
    colSets.push_back(typename Node<T>::SetPtrType(new typename std::set<T>()));
    rowSizeVec.push_back(0);
    colSizeVec.push_back(0);
  }
  for(unsigned int j = 0; j < sqSize; ++j){
    std::vector<typename Node<T>::SetPtrType> tempSets;
    std::vector<unsigned int> tempVec;
    for(unsigned int k = 0; k < sqSize; ++k){
      tempSets.push_back(typename Node<T>::SetPtrType(new typename std::set<T>()));
      tempVec.push_back(0);
    }
    squareSets.push_back(tempSets);
    squareSizeVec.push_back(tempVec);
  }
  
  /********* Fist update the sets for rows, columns and squares ********/
  for(unsigned int row = 0; row < N; ++row){
    for(unsigned int col = 0; col < N; ++col){
      rowSets[row]->insert((*NodeTable)[row][col].getValue());
      colSets[col]->insert((*NodeTable)[row][col].getValue());
      unsigned int sqX = static_cast<unsigned int>(row/sqSize);
      unsigned int sqY = static_cast<unsigned int>(col/sqSize);
      squareSets[sqX][sqY]->insert((*NodeTable)[row][col].getValue());
      if(!(*NodeTable)[row][col].isEmpty()){
	rowSizeVec[row] += 1;
	colSizeVec[col] += 1;
	squareSizeVec[sqX][sqY] += 1;
      }
    }
    if(rowSets[row]->find(zero_value) != rowSets[row]->end()){
      rowSets[row]->erase(rowSets[row]->find(zero_value));
    }
  }
  for(unsigned int col = 0; col < N; ++col){
    if(colSets[col]->find(zero_value) != colSets[col]->end()){
      colSets[col]->erase(colSets[col]->find(zero_value));
    }
  }
  for(unsigned int sqX = 0; sqX < sqSize; ++sqX){
    for(unsigned int sqY = 0; sqY < sqSize; ++sqY){
      if(squareSets[sqX][sqY]->find(zero_value) != squareSets[sqX][sqY]->end()){
	squareSets[sqX][sqY]->erase(squareSets[sqX][sqY]->find(zero_value));
      }
    }
  }
  /********** Now update the set pointers for the individual nodes **********/
  for(unsigned int row = 0; row < N; ++row){
    for(unsigned int col = 0; col < N; ++col){
      unsigned int sqX = static_cast<unsigned int>(row/sqSize);
      unsigned int sqY = static_cast<unsigned int>(col/sqSize);

      (*NodeTable)[row][col].getSquareSet() = squareSets[sqX][sqY];
      (*NodeTable)[row][col].getRowSet()    = rowSets[row];
      (*NodeTable)[row][col].getColSet()    = colSets[col];
      (*NodeTable)[row][col].updatePossibleValues();
      
      /*Check that there are no repeated values in the rows, columns and squares. */
      valid_table = valid_table && ((rowSets[row]->size() == rowSizeVec[row]) && 
				    (colSets[col]->size() == colSizeVec[col]) && 
				    (squareSets[sqX][sqY]->size() == squareSizeVec[sqX][sqY]));
     }
  }
}

template<typename T> std::string Table<T>::tableToString(void){
  std::string out_string(""); 
  for(unsigned int row = 0; row < N; ++row){
    std::string tmp_line("");
    for(unsigned int col = 0; col < N; ++col){
      tmp_line += (boost::lexical_cast<std::string>((*NodeTable)[row][col].getValue())+ " ");
    }
    out_string += (tmp_line + "\n");
  }
  return(out_string);
}

template<typename T>  bool Table<T>::checkIfSolved(std::vector< std::vector<Node<T> > >& recNodeTable){
  bool solved = true;
  for(unsigned int row = 0; row < N; ++row){
    for(unsigned int col = 0; col < N; ++col){
      solved = solved && (!recNodeTable[row][col].isEmpty());
    }
  }
  return(solved);
}

template<typename T> void Table<T>::reverseSinglePossValues(std::vector<std::pair<unsigned int, unsigned int> >& singleValuesSolved){
  for(std::vector<std::pair<unsigned int, unsigned int> >::iterator it = singleValuesSolved.begin(); it != singleValuesSolved.end(); ++it){
    (*NodeTable)[it->first][it->second].setValue(zero_value);
  }
}

template<typename T> bool Table<T>::solveSinglePossValues(std::vector<std::pair<unsigned int, unsigned int> >& singleValuesSolved){
  bool solvable_cells = true;
  bool multi_value_cells = true;
  singleValuesSolved.clear();
  while(solvable_cells && multi_value_cells){
    solvable_cells = false;
    multi_value_cells = false;
    for(unsigned int row = 0; row < N; ++row){
      for(unsigned int col = 0; col < N; ++col){
	if((*NodeTable)[row][col].getPossibleValues()->size() == 1 && (*NodeTable)[row][col].isEmpty()){
	  /* In case there is only one possibility*/
	  solvable_cells = true;
	  T temp_val = static_cast<T>(*((*NodeTable)[row][col].getPossibleValues()->begin()));
	  (*NodeTable)[row][col].setValue(temp_val);
	  singleValuesSolved.push_back(std::pair<unsigned int, unsigned int>(row, col));
	}
	else{
	  if((*NodeTable)[row][col].getPossibleValues()->size() > 1 && (*NodeTable)[row][col].isEmpty()){
	    multi_value_cells = true;
	  }
	}
      }
    }
    updateSets();
  }
  return(multi_value_cells);
}

template<typename T> void Table<T>::recursiveGenValueComb(std::vector<std::pair<unsigned int, unsigned int> > nodeVector, 
							  std::vector<std::pair<unsigned int, unsigned int> > visitedNodeVector, 
							  std::vector<T>& valueVector,
							  std::vector<std::vector<T> >& valueCombinations){
  if(0 == nodeVector.size()){
    typename Node<T>::SetPtrType::element_type local_set(valueVector.begin(), valueVector.end());
    /*If there are repeated values in the row, discard the combination as it gives an invalid table.*/
    if(local_set.size() == valueVector.size()){
      valueCombinations.push_back(valueVector);
    }
    else{
    }
    return;
  }
  std::pair<unsigned int, unsigned int> local_pair(*(nodeVector.rbegin()));
  typename Node<T>::SetPtrType::element_type local_allowed_set((*NodeTable)[local_pair.first][local_pair.second].getPossibleValues()->begin(), 
							       (*NodeTable)[local_pair.first][local_pair.second].getPossibleValues()->end());
  for(typename std::set<T>::iterator it = local_allowed_set.begin(); it != local_allowed_set.end(); ++it){
    typename std::vector<T>::iterator valVecBegin = valueVector.begin();
    valueVector.insert(valVecBegin, *it);
    visitedNodeVector.push_back(local_pair);
    nodeVector.pop_back();
    /////////////////////////////////////////////////////////////////
    recursiveGenValueComb(nodeVector, visitedNodeVector, valueVector, valueCombinations);
    /////////////////////////////////////////////////////////////////
    valVecBegin = valueVector.begin();
    valueVector.erase(valVecBegin);
    visitedNodeVector.pop_back();
    nodeVector.push_back(local_pair);
  }
}

template<typename T>  void Table<T>::recursiveSolve(unsigned int depth){
  ++depth;
  //std::cout << " depth = " + boost::lexical_cast<std::string>(depth) << std::endl;
  if(solution_found){
    std::cout << "Solution found, returning at depth " + boost::lexical_cast<std::string>(depth) << std::endl;
    std::cout << tableToString() << std::endl;
    return;
  }
  
  std::vector<std::pair<unsigned int, unsigned int> > nodeVector; 
  std::vector<std::pair<unsigned int, unsigned int> > visitedNodeVector;
  std::vector<T> valueVector;
  std::vector<std::vector<T> > valueCombinations;
  std::vector<std::pair<unsigned int, unsigned int> > singleValuesSolved;
  unsigned int loop_count = 0;

  /* First find the row with the largest number of already filled squares.*/
  unsigned int c_row          = 0;
  unsigned int max_set_length = 0;
  unsigned int max_set_row    = 0;
  for(typename std::vector< typename Node<T>::SetPtrType >::iterator row_it = rowSets.begin(); row_it != rowSets.end(); ++row_it){
    if((*row_it)->size() == N){
      /* Row is full, continue */
      ++c_row;
      continue;
    }
    else{
      if(max_set_length < (*row_it)->size()){
	max_set_length = (*row_it)->size();
	max_set_row = c_row;
	++c_row;
      }
      else{
	++c_row;
      }
    }
  }
  
  for(unsigned int c_col = 0; c_col < N; ++c_col){
    if((*NodeTable)[max_set_row][c_col].isEmpty()){
      nodeVector.push_back(std::pair<unsigned int, unsigned int>(max_set_row, c_col));
    }
  }
  recursiveGenValueComb(nodeVector, visitedNodeVector, valueVector, valueCombinations);
  for(typename std::vector<std::vector<T> >::iterator v_it = valueCombinations.begin(); v_it != valueCombinations.end(); ++v_it){
    ++loop_count;
    for(unsigned int node_no = 0; node_no < nodeVector.size(); ++node_no){
      T temp_val = static_cast<T>((*v_it)[node_no]);
      (*NodeTable)[nodeVector[node_no].first][nodeVector[node_no].second].setValue(temp_val);
    }      
    updateSets();
    solveSinglePossValues(singleValuesSolved);    
    if(!valid_table){
      for(unsigned int node_no = 0; node_no < nodeVector.size(); ++node_no){
	(*NodeTable)[nodeVector[node_no].first][nodeVector[node_no].second].setValue(zero_value);
      }
      reverseSinglePossValues(singleValuesSolved);
      updateSets();
      continue;
    }
    if(depth < max_depth){
      solution_found = checkIfSolved(*NodeTable);	  
      recursiveSolve(depth);
    }
    if( !solution_found ){
      for(unsigned int node_no = 0; node_no < nodeVector.size(); ++node_no){
	(*NodeTable)[nodeVector[node_no].first][nodeVector[node_no].second].setValue(zero_value);
      }
      reverseSinglePossValues(singleValuesSolved);
      updateSets();
    }
  }
}

template<typename T> void Table<T>::solvePuzzle(void){
  std::vector<std::pair<unsigned int, unsigned int> > singleValuesSolved;
  bool multi_value_cells = solveSinglePossValues(singleValuesSolved);
  if(multi_value_cells){
    unsigned int depth = 0;
    updateSets();
    recursiveSolve(depth);
  }
  else{
    std::cout << tableToString() << std::endl;
  }
}

template<typename T> Table<T>* makeTableFromFile(std::string file_name){
  std::cout << "reading file" << std::endl;
  std::ifstream fin;
  std::string c_line;
  std::vector<std::string> numberRows;
  int line = 0;
  fin.open(file_name.c_str());
  if (fin.is_open()){
    while( fin ){
      getline(fin, c_line);
      if(c_line[0] == '='){
	break;
      }
      numberRows.push_back(c_line);
      ++line;
    }
  }
  fin.close();
  Table<T>* Tbl = new Table<T>(numberRows.size());
  Tbl->updateNodeTable(numberRows);
  Tbl->updateSets();
  return(Tbl);
}

int main(int argc, char* argv[]){
  if(argc < 3){
    std::cerr << "usane ./program_name <file_in> [str, int]" << std::endl;
    exit(-1);
  }
  else{
    if(std::string(argv[2]) == "int"){
      Table<unsigned int>* Tbl = makeTableFromFile<unsigned int>(std::string(argv[1]));
      Tbl->solvePuzzle();
      delete Tbl;
    }
    else{
      if(std::string(argv[2]) == "str"){
	Table<std::string>* Tbl = makeTableFromFile<std::string>(std::string(argv[1]));
	Tbl->solvePuzzle();
	delete Tbl;
      }
      else{
	std::cerr << "Second argument must indicate the type of characters in the puzzle, either int or str" << std::endl;
      }
    }
  }
  return(0);  
}
