#include "basic/pregel-dev.h"
#include <math.h>
#include <set>
#include <limits>

using namespace std;

//input line format: vertexID \t x y
//output line format: 
//	line 1: grid_x grid_y
//	line n>=1: cellID \t point1_x point1_y point2_x point2_y point3_x point3_y...

// Global variables
double EPS;

#define DIMENTION 2
#define ONE_D_SPACE GRID_D1
#define TWO_D_SPACE GRID_D1 * GRID_D2

#define BUF_SIZE 4096
#define BUF_SAFETY_THRESHOLD 4000

int D1_MAX = std::numeric_limits<int>::min();
int D1_MIN = std::numeric_limits<int>::max();
int D2_MAX = std::numeric_limits<int>::min();
int D2_MIN = std::numeric_limits<int>::max();


int GRID_D1;
int GRID_D2;


bool bPrintGrid = true;

int TOTALPTS = 0;
int SEND_COUNT = 0;

typedef struct 
{
	int id1;
	int id2;

} Point;

ibinstream & operator<<(ibinstream & m, const Point & strPoint)
{
	m << strPoint.id1;
	m << strPoint.id2;

	return m;
}

obinstream & operator>>(obinstream & m, Point & strPoint)
{
	m >> strPoint.id1;
	m >> strPoint.id2;

	return m;
}
//==============================================================

typedef struct 
{
	bool bIsCell;
	Point stPointCdnt;
	vector<Point> vctPoints;
} vtxValue;

ibinstream & operator<<(ibinstream & m, const vtxValue & stValue)
{
	m << stValue.bIsCell;
	m << stValue.stPointCdnt;
	m << stValue.vctPoints;
	return m;
}

obinstream & operator>>(obinstream & m, vtxValue & stValue)
{
	m >> stValue.bIsCell;
	m >> stValue.stPointCdnt;
	m >> stValue.vctPoints;
	return m;
}
//==============================================================

class DataVertex: public Vertex<VertexID, vtxValue, Point>
{
private:
	int hashCellId(const Point pnt)
	{
		double dbUnit = EPS / sqrt(DIMENTION);
		GRID_D1 = floor((D1_MAX - D1_MIN) / dbUnit) + 1;
		GRID_D2 = floor((D2_MAX - D2_MIN) / dbUnit) + 1;

		int iNumGrid_d1 = floor((pnt.id1 - D1_MIN) / dbUnit) + 1;
		int iNumGrid_d2 = floor((pnt.id2 - D2_MIN) / dbUnit) + 1;

		//return ((iNumGridY - 1) * GRID_X) + iNumGridX;
		return (iNumGrid_d2 - 1) * ONE_D_SPACE + iNumGrid_d1;
	}

	void sendCellID(const int iCellID)
	{
		Point stCellId;
		stCellId.id1 = iCellID;
		// int iReveiverId = (SEND_COUNT++) % TOTALPTS + 1; 
		send_message(1, stCellId);
	}

public:
	virtual void compute(MessageContainer & messages)
	{
		// Each point get its cell id and send it to other vertices
		if (step_num() == 1) 
		{

			int iCellID = hashCellId(value().stPointCdnt);
			sendCellID(iCellID);

			vote_to_halt();
		}

		// All vertices creates all the cell vertices
		if (step_num() == 2) 
		{
			set<int> setCellIds;
			set<int>::iterator it;
			for (int ii =0 ; ii < messages.size(); ii++)
			{
				setCellIds.insert(messages[ii].id1);
			}

			for (it = setCellIds.begin(); it != setCellIds.end(); it++)
			{
				int iCellVertexId = *it + TOTALPTS;
				// Create cell vertex
				DataVertex *pv = new DataVertex;
				pv->id = iCellVertexId;
				pv->value().bIsCell = true;
				add_vertex(pv);
			}

			wakeAll();

		}

		// Each point get its cell id and send it to the corresponding cell
		if (step_num() == 3)
		{
			if (!value().bIsCell)
			{
				int iCellVertexId = hashCellId(value().stPointCdnt) + TOTALPTS;
				send_message(iCellVertexId, value().stPointCdnt);
				vote_to_halt();
			}

		}

		// Each cell receives it points and save them in the container
		if (step_num() == 4)
		{
			for (int ii = 0; ii < messages.size(); ii++)
			{
				value().vctPoints.push_back(messages[ii]);
			}
			vote_to_halt();
		}

	}
};
//==============================================================

class DataWorker: public Worker<DataVertex>
{
	char buf[BUF_SIZE];	// CAUTION! Need to ensure that the buffer will not overflow, because there may be a lot of points in a cell.

public:
	virtual DataVertex *toVertex(char *pline)
	{
		char *pch;
		pch = strtok(pline, "\t");	// No point id
		DataVertex *pv = new DataVertex;
		// pv->id = atoi(pch);
		pv->id = ++TOTALPTS;
		pv->value().bIsCell = false;

		pch = strtok(NULL, " ");
		int id1 = atoi(pch);
		if (id1 > D1_MAX) 
		{
			D1_MAX = id1;
		}
		if (id1 < D1_MIN)
		{
			D1_MIN = id1;
		}		
		pv->value().stPointCdnt.id1 = id1;

		pch = strtok(NULL, " ");
		int id2 = atoi(pch);
		if (id2 > D2_MAX)
		{
			D2_MAX = id2;
		}
		if (id2 < D2_MIN)
		{
			D2_MIN = id2;
		}
		pv->value().stPointCdnt.id2 = id2;

		return pv;
	}

	virtual void toline(DataVertex *pv, BufferedWriter & writer)
	{
		if (pv->value().bIsCell)
		{
			char *pchBufEnd = buf;
			// Print the width and height of the grid to line 1.
			if (bPrintGrid)
			{
				bPrintGrid = false;
				pchBufEnd += sprintf(pchBufEnd, "$\t%d %d %d %d\n", D1_MIN, D2_MIN, GRID_D1, GRID_D2);
			}

			int iCellID = pv->id - TOTALPTS;
			pchBufEnd += sprintf(pchBufEnd, "%d\t", iCellID);
			for (int ii = 0; ii < pv->value().vctPoints.size(); ii++)
			{
				pchBufEnd += sprintf(pchBufEnd, "%d %d ", pv->value().vctPoints[ii].id1, pv->value().vctPoints[ii].id2);
				if (pchBufEnd - buf > BUF_SAFETY_THRESHOLD)
				{
					writer.write(buf);
					memset(buf, 0, BUF_SIZE);
					pchBufEnd = buf;
				}			
			}
			pchBufEnd += sprintf(pchBufEnd, "\n");
			writer.write(buf);
		}
	}

};
//=============================================================

void pregel_pre_data(string in_path, string out_path, double dbEPS)
{
	::EPS = dbEPS;

	WorkerParams param;
	param.input_path = in_path;
	param.output_path = out_path;
	param.force_write = true;
	param.native_dispatcher = false;
	DataWorker worker;

	worker.run(param);
}
