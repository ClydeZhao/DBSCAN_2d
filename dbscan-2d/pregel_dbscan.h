#include "basic/pregel-dev.h"
#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;

// Input format: 
//	line 1: GRID_D1 GRID_D2
//	line n>=1: cellID \t point1_x point1_y point2_x point2_y point3_x point3_y...
// Output format:
// clusterID \t x y

// Global variables
#define DIMENTION 2

#define BUF_SIZE 4096
#define BUF_THRSHOLD 4000

double EPS;
double p;
int MinPts;

int D1_MIN;
int D2_MIN;


int GRID_D1;
int GRID_D2;



enum MessageType
{
	REQ,
	CELL,
	MIN,
	POINT
};

ibinstream & operator<<(ibinstream & m, const MessageType & eType)
{
	m << eType;
	return m;
}

obinstream & operator>>(obinstream & m, MessageType & eType)
{
	m >> eType;
	return m;
}
//=========================================================================

struct Point
{
	int id1;
	int id2;
	bool bIsCore;
	bool bIsEarsed;
};

ibinstream & operator<<(ibinstream & m, const Point & stPoint)
{
	m << stPoint.id1;
	m << stPoint.id2;
	m << stPoint.bIsCore;
	m << stPoint.bIsEarsed;
	return m;
}

obinstream & operator>>(obinstream & m, Point & stPoint)
{
	m >> stPoint.id1;
	m >> stPoint.id2;
	m >> stPoint.bIsCore;
	m >> stPoint.bIsEarsed;
	return m;
}
//=========================================================================

struct Cell
{
	double dbd1;
	double dbd2;
	double dbLength;
	int iClusterId;
	bool bIsCore;
	vector<int> vtAdjList;
	vector<Point> vtDataPoints;	
};

ibinstream & operator<<(ibinstream & m, const Cell & stCell)
{
	m << stCell.dbd1;
	m << stCell.dbd2;
	m << stCell.dbLength;
	m << stCell.iClusterId;
	m << stCell.bIsCore;
	m << stCell.vtAdjList;
	m << stCell.vtDataPoints;
	return m;
}

obinstream & operator>>(obinstream & m, Cell & stCell)
{
	m >> stCell.dbd1;
	m >> stCell.dbd2;
	m >> stCell.dbLength;
	m >> stCell.iClusterId;
	m >> stCell.bIsCore;
	m >> stCell.vtAdjList;
	m >> stCell.vtDataPoints;
	return m;
}
//=======================================================================

struct Message
{
	MessageType type;
	int iSenderId;
	Cell stCell;
	Point stPoint;
};

ibinstream & operator<<(ibinstream & m, const Message & stM)
{
	m << stM.type;
	m << stM.iSenderId;
	m << stM.stCell;
	m << stM.stPoint;
	return m;
}

obinstream & operator>>(obinstream & m, Message & stM)
{
	m >> stM.type;
	m >> stM.iSenderId;
	m >> stM.stCell;
	m >> stM.stPoint;
	return m;
}
//=======================================================================

class TreeNode
{
private:
	double dbd1;
	double dbd2;
	double dbLength;
	int iCount;

	int hash_data_point(const Point stPnt, const Cell stCell)
	{
		double dbSubCellLength = stCell.dbLength / 2;
		int iSeq = 1;
		if (stPnt.id1 > stCell.dbd1 + dbSubCellLength)
		{
			iSeq += 1;
		}
		if (stPnt.id2 > stCell.dbd2 + dbSubCellLength)
		{
			iSeq += 2;
		}
		return iSeq;

	}

	vector<Cell> split_cell(const Cell stCell, const bool bCoreTree)
	{
		vector<Cell> vtSubCells;
		double dbSubCellLength = stCell.dbLength / 2;
		vector<Point> vtDataPoints = stCell.vtDataPoints;
		int iNumOfSubCells = pow(2, DIMENTION);

		for (int ii = 1; ii <= iNumOfSubCells; ii++)
		{
			Cell stSubCell;
			stSubCell.dbLength = dbSubCellLength;
			if (ii % 2 == 1)
			{
				stSubCell.dbd1 = stCell.dbd1;
			}
			else
			{
				stSubCell.dbd1 = stCell.dbd1 + dbSubCellLength;
			}

			if (ii % 4 <= 2)
			{
				stSubCell.dbd2 = stCell.dbd2;
			}
			else
			{
				stSubCell.dbd2 = stCell.dbd2 + dbSubCellLength;
			}
			vtSubCells.push_back(stSubCell);
		}
		

		// Hash data points to sub cells
		for (int ii = 0; ii < vtDataPoints.size(); ii++)
		{
			if (bCoreTree)
			{
				if (vtDataPoints[ii].bIsCore)
				{
					int subCellNum = hash_data_point(vtDataPoints[ii], stCell);
					vtSubCells[subCellNum-1].vtDataPoints.push_back(vtDataPoints[ii]);
				}
			}
			else
			{
				int subCellNum = hash_data_point(vtDataPoints[ii], stCell);
				vtSubCells[subCellNum-1].vtDataPoints.push_back(vtDataPoints[ii]);
			}
		}

		return vtSubCells;
	}


public:
	vector<TreeNode *> vtpChildren;

	TreeNode();

	TreeNode(const Cell stCell)
	{
		this->dbd1 = stCell.dbd1;
		this->dbd2 = stCell.dbd2;
		this->dbLength = stCell.dbLength;
		iCount = stCell.vtDataPoints.size();
	}

	// ~TreeNode();
	// {
	// 	for (int ii = 0; ii < vtpChildren.size(); ii++)
	// 	{
	// 		delete vtpChildren[ii];
	// 	}
	// }

	void Build_tree(TreeNode *pParent, const Cell stCell, bool bCoreTree = false)
	{
		if (stCell.dbLength < (2 * EPS * p / sqrt(DIMENTION) ))
		{
			return;
		}
		else
		{
			vector<Cell> vtSubCells = split_cell(stCell, bCoreTree);
			int iNumOfSubCells = pow(2, DIMENTION);
			for (int ii = 0; ii < iNumOfSubCells; ii++)
			{
				if (vtSubCells[ii].vtDataPoints.size() != 0)
				{
					TreeNode *pChild = new TreeNode(vtSubCells[ii]);
					pParent->vtpChildren.push_back(pChild);
					Build_tree(pChild, vtSubCells[ii], bCoreTree);
				}
			}
		}
	}

	vector<double *> Get_corners() const
	{
		vector<double *> vcCorners;

		int inumOfCorners = pow(2, DIMENTION);
		for (int ii = 1; ii <= inumOfCorners; ii++)
		{
			double * dbCorner = new double[DIMENTION];
			if (ii % 2 == 1)
			{
				dbCorner[0] = this->dbd1;

			}
			else
			{
				dbCorner[0] = this->dbd1 + this->dbLength;

			}

			if (ii <= 2)
			{
				dbCorner[1] = this->dbd2;
			}
			else
			{
				dbCorner[1] = this->dbd2 + this->dbLength;
			}

			vcCorners.push_back(dbCorner);
		}

		return vcCorners;
	}

	int Get_count() const 
	{
		return this->iCount;
	}

	bool Is_leaf() const
	{
		if (this->vtpChildren.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	// For debug
	void Print_corners() const	
	{
	 	std::vector<double *> vtCorners = this->Get_corners();
		printf("%f\t%f\t%f\t%d\n", this->dbd1, this->dbd2, this->dbLength, this->iCount);
		printf("{(%f, %f) (%f, %f) (%f, %f) (%f, %f)}\n", vtCorners[0][0], vtCorners[0][1], vtCorners[1][0], vtCorners[1][1], vtCorners[2][0], vtCorners[2][1], vtCorners[3][0], vtCorners[3][1]);
		for (int ii = 0; ii < vtCorners.size(); ii++)
		{
			delete[] vtCorners[ii];
		}
	}

};
//==========================================================================

class CellVertex: public Vertex<VertexID, Cell, Message>
{
private:
	void label_core_cell()
	{
		value().bIsCore = true;
		for (int ii = 0; ii < value().vtDataPoints.size(); ii++)
		{
			value().vtDataPoints[ii].bIsCore = true;
		}
	}


	int * get_dim_array(const int id, const int dimSpace, const int range)
	{
		int * dimArr = new int[5];
		dimArr[0] = -2 * dimSpace;
		dimArr[1] = -1 * dimSpace;
		dimArr[2] = 0;
		dimArr[3] = dimSpace;
		dimArr[4] = 2 * dimSpace;

		if (ceil(double(id) / dimSpace) == 1)	
		{
			dimArr[0] = dimArr[1] = 0;
		}
		if (ceil(double(id) / dimSpace) == 2)
		{
			dimArr[0] = 0;
		}
		if (ceil(double(id) / dimSpace) == range)
		{
			dimArr[3] = dimArr[4] = 0;
		}
		if (ceil(double(id) / dimSpace) == range - 1)
		{
			dimArr[0] = 0;
		}
		return dimArr;
	}

	vector<int> get_eps_neighor_id(const int center)
	{
		int index = 0;
		vector<int> vtResult;
		int OneDSpace = GRID_D1;
		int TwoDSpace = OneDSpace * GRID_D2;

		int TwoDId = center;

		int OneDId = TwoDId	% OneDSpace;
		if (OneDId == 0)
		{
			OneDId = OneDSpace;
		}

		int * TwoDArr = get_dim_array(TwoDId, OneDSpace, GRID_D2);
		int * OneDArr = get_dim_array(OneDId, 1, GRID_D1);

		for (int ii = 0; ii < 5; ii++)
		{
			if (OneDArr[ii] == 0 && ii != 2)
			{
				continue;
			}
			for (int jj = 0; jj < 5; jj++)
			{
				if (TwoDArr[jj] == 0 && jj != 2)
				{
					continue;
				}
				index = center + OneDArr[ii] + TwoDArr[jj];
				if (index > 0 && index <= TwoDSpace)
				{
					vtResult.push_back(index);
				}									
			}
		}

		delete[] TwoDArr;
		delete[] OneDArr;

		return vtResult;
	}

	void send_req_to_eps_neighbor()
	{
		Message stReqMsg;
		stReqMsg.type = REQ;
		stReqMsg.iSenderId = id;
		vector<int> vtEpsNeighborIds = get_eps_neighor_id(id);
		for (int ii = 0; ii <  vtEpsNeighborIds.size(); ii++)
		{
			send_message(vtEpsNeighborIds[ii], stReqMsg);

		}

	}

	double compute_distance(const Point stPnt, const double * corner)
	{
		// return sqrt(pow(stPnt.ix - corner.first, 2) + pow(stPnt.iy - corner.second, 2));
		return sqrt(pow(stPnt.id1 - corner[0], 2) + pow(stPnt.id2 - corner[1], 2));
	}

	void delete_corner_vector(vector<double *> &v)
	{
		for (int ii = 0; ii < v.size(); ii++)
		{
			delete[] v[ii];
		}
	}

	bool is_disjointed_with_ball(const Point stPnt, const TreeNode *pRoot)
	{
		vector<double  *> vtCorners = pRoot->Get_corners();
		int inumOfCorners = pow(2, DIMENTION);

		for (int ii = 0; ii < inumOfCorners; ii++)
		{
			if (compute_distance(stPnt, vtCorners[ii]) < EPS)
			{
				return false;
			}
		}

		delete_corner_vector(vtCorners);
		return true;
	}

	bool is_fully_covered_by_p_ball(const Point stPnt, const TreeNode *pRoot)
	{
		double dbDiameter = (1 + p) * EPS;
		vector<double  *> vtCorners = pRoot->Get_corners();
		int inumOfCorners = pow(2, DIMENTION);

		for (int ii = 0; ii < inumOfCorners; ii++)
		{
			if (compute_distance(stPnt, vtCorners[ii]) > dbDiameter)
			{
				return false;
			}
		}

		delete_corner_vector(vtCorners);
		return true;
	}

	void recursive_count(const Point stPnt, const TreeNode *pRoot, int & ans)
	{
		if (is_disjointed_with_ball(stPnt, pRoot))
		{
			//ans += 0;
			return;
		}
		if (is_fully_covered_by_p_ball(stPnt, pRoot))
		{
			ans += pRoot->Get_count();
			return;
		}
		if (pRoot->Is_leaf())
		{
			ans += pRoot->Get_count();
			return;
		}
		else
		{
			for (int ii = 0; ii < pRoot->vtpChildren.size(); ii++)
			{
				recursive_count(stPnt, pRoot->vtpChildren[ii], ans);
			}
		}
	}

	int approximate_range_count(const Point stPnt, const TreeNode *pRoot)
	{
		int ans = 0;

		recursive_count(stPnt, pRoot, ans);

		return ans;

	}

	void hashmin_broadcast(const int min)
	{
		for (int ii = 0; ii < value().vtAdjList.size(); ii++)
		{
			Message stMinMsg;
			stMinMsg.type = MIN;
			stMinMsg.iSenderId = min;

			send_message(value().vtAdjList[ii], stMinMsg);
		}
	}

	//debug
	// void print_point(const Point p) const
	// {
	// 	printf("(%d, %d)  ", p.ix, p.iy);
	// }

	// void print_cell_info() const
	// {
	// 	printf("CellID: %d, X: %f  Y: %f\n", id, value().dbx, value().dby);
	// 	for (int ii = 0; ii < value().vtDataPoints.size(); ii++)
	// 	{
	// 		print_point(value().vtDataPoints[ii]);
	// 	}
	// 	printf("\n");
	// }

public:
	virtual void compute(MessageContainer & messages)
	{
		// Step 1: If |Cell| > MinPts, label cell as core; else, send request to Epsilon-neighbors.
		if (step_num() == 1)
		{
			if (id == 0)
			{
				vote_to_halt();
				return;
			}
			if (value().vtDataPoints.size() >= MinPts)
			{
				// value().bIsCore = true;
				label_core_cell();
			}
			else
			{
				send_req_to_eps_neighbor();
			}
			vote_to_halt();
		}

		// Step 2: On receiving request message, send cell information to requesters
		if (step_num() == 2)
		{
			Message stInfoMsg;
			stInfoMsg.type = CELL;
			stInfoMsg.iSenderId = id;
			stInfoMsg.stCell = value();

			for (int ii = 0; ii < messages.size(); ii++)
			{
				if (messages[ii].type == REQ)
				{
					send_message(messages[ii].iSenderId, stInfoMsg);
				}
			}
			vote_to_halt();
		}

		// Step 3: On receiving cell info messages, non-core cell contruct tree and do approximate range count for labelling.
		if (step_num() == 3)
		{

			vector<TreeNode> vtNeighborTrees;
			for (int ii = 0; ii < messages.size(); ii++)
			{
				if (messages[ii].type == CELL)
				{
					TreeNode root(messages[ii].stCell);
					root.Build_tree(&root, messages[ii].stCell);
					vtNeighborTrees.push_back(root);

				}

			}
						
			int iCount;
			int iTotal;
			for (int ii = 0; ii < value().vtDataPoints.size(); ii++)
			{
				iTotal = 1;
				for (int jj = 0; jj < vtNeighborTrees.size(); jj++)
				{
					iCount = approximate_range_count(value().vtDataPoints[ii], &(vtNeighborTrees[jj]) );
					iTotal += iCount;	
		
				}

				if (iTotal >= MinPts)
				{
					value().vtDataPoints[ii].bIsCore = true;
					value().bIsCore = true;
					break;
				}
			}
			wakeAll();
		}

		// Build a graph connecting all the core cells
		// Step 4: Each core cell request its epsilon neighbors' information, so as to do approximate count in the next  step
		if (step_num() == 4)
		{
			if (value().bIsCore)
			{
				send_req_to_eps_neighbor();
			}
			else
			{
				vote_to_halt();
			}
		}
		// Step 5: Respond to request.
		if (step_num() == 5)
		{
			if (value().bIsCore)
			{
				Message stInfoMsg;
				stInfoMsg.type = CELL;
				stInfoMsg.iSenderId = id;
				stInfoMsg.stCell = value();

				for (int ii = 0; ii < messages.size(); ii++)
				{
					if (messages[ii].type == REQ)
					{
						send_message(messages[ii].iSenderId, stInfoMsg);
					}
				}
				vote_to_halt();
			}
		}

		// Step 6: Core cells receive epsilon neighbors' information, and do approximate count. If count is positive, add the neighbor to adjacency lists
		if (step_num() == 6)
		{
			int iCount = 0;
			for (int ii = 0; ii < messages.size(); ii++)
			{
				if (messages[ii].type == CELL && messages[ii].stCell.bIsCore == true)
				{
					TreeNode root(messages[ii].stCell);
					root.Build_tree(&root, messages[ii].stCell, true);	// Core tree
					for (int jj = 0; jj < value().vtDataPoints.size(); jj++)
					{
						if (value().vtDataPoints[jj].bIsCore)
						{
							iCount = approximate_range_count(value().vtDataPoints[jj], &root);

							if (iCount > 0)
							{
								value().vtAdjList.push_back(messages[ii].iSenderId);
								break;
							}
						}
					}
				}
			}
			wakeAll();
		}

		// Step 7 - 10: Send boarder points to a near core cell
		if (step_num() == 7)
		{
			// Non-core cells
			if (value().bIsCore == false)
			{
				send_req_to_eps_neighbor();
			}
			vote_to_halt();

		}

		if (step_num() == 8)
		{
			if (value().bIsCore)
			{
				Message stInfoMsg;
				stInfoMsg.type = CELL;
				stInfoMsg.iSenderId = id;
				stInfoMsg.stCell = value();

				for (int ii = 0; ii < messages.size(); ii++)
				{
					if (messages[ii].type == REQ)
					{
						send_message(messages[ii].iSenderId, stInfoMsg);
					}
				}
			}
			vote_to_halt();
		}

		if (step_num() == 9)
		{
			if (value().bIsCore == false) // make sure that this is a non-core cell
			{
				vector<Point> &vtDataPoints = value().vtDataPoints;
				vector<pair<int, TreeNode> > vtNeighborTrees;
				for (int ii = 0; ii < messages.size(); ii++)
				{
					if (messages[ii].type == CELL)
					{
						TreeNode root(messages[ii].stCell);
						root.Build_tree(&root, messages[ii].stCell, true);	// Core tree
						vtNeighborTrees.push_back(make_pair(messages[ii].iSenderId, root));
					}
				}

				int iCount = 0;
				vector<int> vtToErease;
				for (int ii = 0; ii < vtDataPoints.size(); ii++)
				{
					for (int jj = 0; jj < vtNeighborTrees.size(); jj++)
					{
						iCount = approximate_range_count(vtDataPoints[ii], &(vtNeighborTrees[jj].second));
						if (id == 292)
						{
							printf("%d %d\n", iCount, vtNeighborTrees[jj].first);
						}
						if (iCount > 0)
						{
							Message stPntMsg;
							stPntMsg.type = POINT;
							stPntMsg.stPoint = vtDataPoints[ii];
							send_message(vtNeighborTrees[jj].first, stPntMsg);
							vtDataPoints[ii].bIsEarsed = true;
							
							break;
						}
					}
				}
			}
			wakeAll();	//In case that there is no boarder point
		}

		if (step_num() == 10)
		{
			if (value().bIsCore)
			{
				for (int ii = 0; ii < messages.size(); ii++)
				{
					if (messages[ii].type == POINT)
					{
						value().vtDataPoints.push_back(messages[ii].stPoint);
					}
				}
			}
			wakeAll();
		}

		// Hashmin: find connected components (clustering)
		if (step_num() == 11)
		{
			if (value().bIsCore)
			{
				int min = id;
				for (int ii = 0; ii < value().vtAdjList.size(); ii++)
				{
					if (min > value().vtAdjList[ii])
					{
						min = value().vtAdjList[ii];
					}
				}

				value().iClusterId = min;
				hashmin_broadcast(min);
				
			}			
			vote_to_halt();
		}

		if (step_num() > 11)
		{
			int min = messages[0].iSenderId;
			for (int ii = 1; ii < messages.size(); ii++)
			{
				if (min > messages[ii].iSenderId)
				{
					min = messages[ii].iSenderId;
				}
			}
			if (min < value().iClusterId)
			{
				value().iClusterId = min;
				hashmin_broadcast(min);
			}
			vote_to_halt();
		}
	}
};

class DbscanWorker: public Worker<CellVertex>
{
	char buf[BUF_SIZE];

public:
	virtual CellVertex *toVertex(char *pline)
	{
		char *pch;
		pch = strtok(pline, "\t");
		// Grid information
		if (*pch == '$')
		{
			//cout << "toVertex： grid" << endl;
			pch = strtok(NULL, " ");
			D1_MIN = atoi(pch);
			pch = strtok(NULL, " ");
			D2_MIN = atoi(pch);

			pch = strtok(NULL, " ");
			GRID_D1 = atoi(pch);
			pch = strtok(NULL, " ");
			GRID_D2 = atoi(pch);

			CellVertex *pv = new CellVertex;
			pv->id = 0;
			pv->value().bIsCore = false;

			return pv;

		}

		else
		{
			//cout << "toVertex： cell" << endl;
			double dbUnit = EPS / sqrt(DIMENTION);

			CellVertex *pv = new CellVertex;
			int grid_id = atoi(pch);
			pv->id = grid_id;
			pv->value().bIsCore = false;
			pv->value().iClusterId = 0;

			int OneDSpace = GRID_D1;
			int TwoDSpace = GRID_D2 * OneDSpace;

			int OneDId = grid_id % OneDSpace;
			if (OneDId == 0)
			{
				OneDId = OneDSpace;
			}

			int iTwoDPos = ceil((double)grid_id / OneDSpace);
			int iOneDPos = OneDId;

			double dbd1 = D1_MIN + (iOneDPos - 1) * dbUnit;
			double dbd2 = D2_MIN + (iTwoDPos - 1) * dbUnit; 

			pv->value().dbd1 = dbd1;
			pv->value().dbd2 = dbd2;
			pv->value().dbLength = dbUnit;


			// CAUTION! Data formate is asumed to be correct. Wrong formate will lead to crash?
			pch = strtok(NULL, " ");
			while(pch != NULL)
			{
				Point stPnt;
				int id1 = atoi(pch);
				pch = strtok(NULL, " ");
				int id2 = atoi(pch);

				stPnt.id1 = id1;
				stPnt.id2 = id2;

				stPnt.bIsCore = false;
				stPnt.bIsEarsed = false;
				pv->value().vtDataPoints.push_back(stPnt);

				pch = strtok(NULL, " ");
			}

			return pv;

		}
	}

	// CAUTION!!! Buffer overflow!!!!!!
	virtual void toline(CellVertex *pv, BufferedWriter & writer)
	{
		*buf = '\0';
		char *pchBufEnd = buf;
		for (int ii = 0; ii < pv->value().vtDataPoints.size(); ii++)
		{
			if (pv->value().vtDataPoints[ii].bIsEarsed == false)
			{
				pchBufEnd += sprintf(pchBufEnd, "%d %d\t%d\n", pv->value().vtDataPoints[ii].id1, pv->value().vtDataPoints[ii].id2,  pv->value().iClusterId);
				if (pchBufEnd - buf > BUF_THRSHOLD)
				{
					writer.write(buf);
					memset(buf, 0, BUF_SIZE);
					pchBufEnd = buf;
				}
			}

		}
		//sprintf(buf, "%d %d\n", pv->id, pv->value().iClusterId);
		writer.write(buf);
	}

};

void pregel_dbscan(string in_path, string out_path, double dbEPS, int iMinPts, double dbp)
{
	EPS = dbEPS;
	MinPts = iMinPts;
	p = dbp;

	WorkerParams param;
	param.input_path = in_path;
	param.output_path = out_path;
	param.force_write = true;
	param.native_dispatcher = false;
	DbscanWorker worker;

	worker.run(param);
}
