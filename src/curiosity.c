#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
/***************************************************************** 
 * a minimalistic mini-game where one has to guide Curiosity     *
 * (robot on mars) through a perilous labyrinth in search of     *
 * intelligent (i.e. feline) life                                *
 *                                                               *  
 * Bart Pieters 2024                                             *            
 *****************************************************************/ 

// some structs
typedef struct { // maze struct with size (Nx, Ny) connections to cells below and right, and whether the cell was visited
	int Nx, Ny; 		// regular grid of cells
	char *down, *right; // connections down/right (0 no connection, 1 connection)
	char *visited; 		// visited (used for maze generation)
} Maze;

typedef struct { // curiosity is in a maze (cell index n) and is oriented (o)
	int n; 		// location (cell index)
	char o; 	// orientation: facing u(p) d(own) r(ight) l(eft)
} curiosity;

// a route is a sequence of visited cells and a bunch of a characters to draw the trail
typedef struct {
	int *n; // cells
	int *pd; // array of route path drawing codes, see below
	int N;
} route;

// route path drawing codes 
#define R_horizontal 	0
#define R_vertical 		1
#define R_right_down 	2
#define R_right_up 		3
#define R_left_down 	4
#define R_left_up 		5
#define R_down_left 	6
#define R_down_right 	7
#define R_up_left 		8
#define R_up_right 		9
#define R_right_left 	10
#define R_left_right 	11
#define R_down_up 		12
#define R_up_down 		13
#define R_no_path 		14


// tables with drawing elements for curiosity and the maze
// for portability better stick with plain ASCII
//#define UNICODEMAZE
//#define UNICODEROBOT
#ifdef UNICODEROBOT
// curiosity drawing vs orientation
char *Cu="🡱"; // up
char *Cd="🡳"; // down
char *Cl="🡰"; // left
char *Cr="🡺"; // right
#else
// curiosity drawing vs orientation
char *Cu="^"; // up
char *Cd="v"; // down
char *Cl="<"; // left
char *Cr=">"; // right
#endif //ifdef UNICODEROBOT
#ifdef UNICODEMAZE
char *paths[]={
	"━",			// horizontal
	"┃",			// vertical
	"┓",			// right-down 
	"┛",			// right-up
	"┏",			// left-down
	"┗",			// left-up
	"┛",			// down-left
	"┗",			// down-right
	"┓",			// up-left
	"┏",			// up-right
	"╸",			// right-left
	"╺",			// left-right
	"╹",			// down-up
	"╻",			// up-down
	" "				// no-path
};

char *innerborders[]={
	"───┼", //0
	"   ├", //1
	"───┬", //2
	"   ┌", //3
	"───┤", //4
	"   │", //5
	"───┐", //6
	"   ╷", //7
	"───┴", //8
	"   └", //9
	"────", //10
	"   ╶", //11
	"───┘", //12
	"   ╵", //13
	"───╴", //14
	"   ",  //15
	"│"     //16
};

char *outerborders[]={
	"┏",    //0 0-3: top border
	"━━━━", //1
	"━━━┯", //2
	"━━━┓", //3
	"┃",    //4 left/right
	"┠",	//5 left w inner
	"───┨", //6 right w inner
	"┗",    //7 lower left corner
	"━━━┷", //8 lower w inner
	"━━━┛", //9 lower right corner
};
#else // do not use unicode
char *paths[]={
	"-",			// horizontal
	"|",			// vertical
	"v",			// right-down 
	"^",			// right-up
	"v",			// left-down
	"^",			// left-up
	"<",			// down-left
	">",			// down-right
	"<",			// up-left
	">",			// up-right
	"-",			// right-left
	"-",			// left-right
	"|",			// down-up
	"|",			// up-down
	" "				// no-path
};

char *innerborders[]={
	"---+", //0
	"   +", //1
	"---+", //2
	"   +", //3
	"---+", //4
	"   |", //5
	"---+", //6
	"   +", //7
	"---+", //8
	"   +", //9
	"----", //10
	"   +", //11
	"---+", //12
	"   +", //13
	"---+", //14
	"   ",  //15
	"|"     //16
};

char *outerborders[]={
	"+",    //0 0-3: top border
	"----", //1
	"---+", //2
	"---+", //3
	"|",    //4 left/right
	"+",	//5 left w inner
	"---+", //6 right w inner
	"+",    //7 lower left corner
	"---+", //8 lower w inner
	"---+", //9 lower right corner
};
#endif // UNICODEMAZE

/* curiosity may:
 * - Aquire the target (success)
 * - Crash in a wall (not so good)
 * - Get lost (also not good)
 * - Receive an illegal instruction and malfunction
 * For each these cases there are several possible messages you may get:
 * The messages may contain a pause at some point. To insert a pause put a tag in a message formatted as "<P12>"
 * where the number (12 in this case) refers to the number of seconds to wait at the insertion point of the tag
 */

int NMcrash=3;
//   ----------------------------------------
char *Mcrash[] = {
	"Oh No! You killed Curiosity\n",
	"Curiosity crashed against a wall.\nthanks a lot, dumbass!\n",
	"You call that programming! You drove\nCuriosity into a wall!"
};
int NMlost=2;
//   ----------------------------------------
char *Mlost[] = {
	"Curiosity is lost, dwelling the martian\nlabyrinth\n",
	"You lost curiosity!\n"
};
int NMsuccess=5;
//   ----------------------------------------
char *Msuccess[] = {
	"! robot finds kitten !\n",
	"False alarm! it was just a rock\n",
	"It is a cat.<P2>\nOh No! Curiosity killed the cat!\n",
	"It's a cat, Schrödinger's cat.\nIt put curiosity in a box!\n",
	"It is life but not intelligent (i.e.\nnot a cat), so never mind\n"
};
int NMillegal=3;
//   ----------------------------------------
char *Millegal[] = {
	"Illegal instruction at address @152637.\nLost connection to curiosity!\n",
	"Segmentation fault.\nLost connection to curiosity!\n",
	"Illegal instruction at address @6353298.\nYou are staring at a blue screen\n"
};
int NMhitnrun=3;
//   ----------------------------------------
char *Mhitnrun[] = {
	"Curiosity drove right over a delicate\nfeline life form.\n",
	"Curiosity to ground control:\nI hit a cat!<P2>\nI think nobody saw it was me.\n",
	"Curiosity to ground control:\nI hit a cat! It looks angry...<P2>\nLost connection to Curiosity!\n"
};


// introduction of the story
char * Intro =          
"----------------------------------------\n"\
" ,--.  ,--.  ,---.   ,---.    ,---.   \n"\
" |  ,'.|  | /  O  \\ '   .-'  /  O  \\  \n"\
" |  |' '  ||  .-.  |`.  `-. |  .-.  | \n"\
" |  | `   ||  | |  |.-'    ||  | |  | \n"\
" `--'  `--'`--' `--'`-----' `--' `--' \n"\
"             NEEDS YOU!               \n"\
"----------------------------------------\n"\
"We, from SETI, have possibly detected an \n"\
"intelligent life form on Mars! This is \n"\
"HUGE! We immediately sent Curiosity to \n"\
"investigate. However, a perilous maze \n"\
"stands between Curiosity and the target.\n"\
"As this mission is of the utmost \n"\
"importance we implore onto you, as the\n"\
"most qualified person for the job, to \n"\
"write the code that brings Curiosity \n"\
"safely to the target.\n"\
"----------------------------------------\n";

// how to play
char * Situation = 
"----------------------------------------\n"\
"Below you will see the situation, with \n"\
"Curiosity marked with %s (note the \n"\
"orientation) and the object marked with \n"\
"O.\n\n"\
"----------------------------------------\n";

char * Guide = 
"----------------------------------------\n"\
"Your Job:\n"\
"Give Curiosity the right instructions \n"\
"to safely arrive at the object. All\n"\
"navigational instructions consist of an\n"\
"optional rotation and one step forwards.\n"\
"The available instructions are:\n"\
"  - \"f\" - forward: no rotation, only \n"\
"          move forwards\n"\
"  - \"b\" - backward: rotate 180 degrees \n"\
"         and move forwards\n"\
"  - \"l\" - left: rotate 90 degrees \n"\
"         counter clockwise and move\n"\
"         forwards.\n"\
"  - \"r\" - right: rotate 90 degrees \n"\
"         clockwise and move forwards\n"\
"You need to provide the complete sequence \n"\
"of commands in one go!\n"\
"----------------------------------------\n";
	
// code	 
int Random(int rmin, int rmax)
{
	int R;
	R=rand();
	return rmin + (int) (1.0*(rmax-rmin+1) * (R-1.0)/RAND_MAX );
}
void InitRandom(void)
{
	time_t curtime;
	
	time( &curtime );
	srand( (unsigned int) curtime );
	Random(0, 1000);
}

// mazes use row major ordering
#define ROW(n,M) ((M->Nx) == 0 ? 0 : ((n)/(M->Nx)))
#define COL(n,M) ((M->Nx) == 0 ? 0 : ((n)%(M->Nx)))
#define INDEX(r,c,M) ((r)*(M->Nx)+(c))


/* we use Wilson's algorithm
 * Meaning we do random walks from some arbitrary cell ultil we hit a
 * cell that was already visited, where we delete loops (i.e. we do not 
 * allow a random walk to visit itself.
 */
// 0: right
// 1: down
// 2: left
// 3: up
// choose a random *and* possible direction 
int RandomDirConstrained(int r, int c, Maze *M)
{
	int d=-1;
	while (d<0)
	{
		d=Random(0, 3);
		if ((c>=M->Nx-1)&&(d==0)) // no right
			d=-1;
		if ((c<=0)&&(d==2)) // no left
			d=-1;
		if ((r>=M->Ny-1)&&(d==1)) // no down
			d=-1;
		if ((r<=0)&&(d==3)) // no up
			d=-1;
	}
	return d;
}	
// delete a loop
int DeleteLoop(int *Path, int i, Maze *M)
{
	int s;
	s=Path[i];
	i--;
	while ((Path[i]!=s)&&(i>0))
	{
		M->visited[Path[i]]=0;
		i--;
	}
	return i+1;
}
// add a path to the maze
void AddPath(int *Path, int i, Maze *M)
{
	int r,c, r2, c2, dr, dc;
		
	M->visited[Path[i]]=1;
	r=ROW(Path[i],M);
	c=COL(Path[i],M);
	i--;
	for (;i>=0;i--)
	{
		M->visited[Path[i]]=1;
		r2=ROW(Path[i],M);
		c2=COL(Path[i],M);
		dr=r2-r;
		dc=c2-c;
		
		if (dr==1) // down
			M->down[Path[i+1]]=1;
		else if (dr==-1) // up
			M->down[Path[i]]=1;
		else if (dc==1) // right
			M->right[Path[i+1]]=1;
		else if (dc==-1) // left
			M->right[Path[i]]=1;
		r=r2;
		c=c2;
	}	
}
// make a random walk to add a path to the maze
void RandomWalk(Maze *M, int n)
{
	int i, d;
	int *Path;
	int r, c;
	
	if (M->visited[n])
		return;		
	M->visited[n]=2;// mark as 2
	
	// not saving memory here, but it should be small anyway
	if ((Path=malloc(M->Nx*M->Ny*sizeof(int)))==NULL)
		return;
		
	r=ROW(n,M);
	c=COL(n,M);
	i=1;
	Path[0]=n;
	while(1)
	{
		d=RandomDirConstrained(r, c, M);
		switch(d)
		{
			case 0:// go right				
				c++;
				break;
			case 1:// go down	
				r++;
				break;
			case 2:// go left	
				c--;
				break;
			default:// go up	
				r--;
				break;
		}
		n=INDEX(r,c,M);
		Path[i]=n;
		if (M->visited[n]==1)
		{
			// add path
			AddPath(Path, i, M);
			free(Path);
			return;
		}
		else if (M->visited[n]==2)
		{
			i=DeleteLoop(Path, i, M);
			r=ROW(n,M);
			c=COL(n,M);
		}
		else
		{
			M->visited[n]=2;
			i++;
		}
	}
}
// select a possible cell to start a random walk
int SelectStart(Maze *M)
{
	int *NV;
	int i, j;
	// not saving memory here, but it should be small anyway
	if ((NV=malloc(M->Nx*M->Ny*sizeof(int)))==NULL)
		return -1;
	j=0;
	for (i=0;i<M->Nx*M->Ny;i++)
		if (M->visited[i]==0)
			NV[j++]=i;
	if (j==0)
	{
		free(NV);
		return -1;
	}
	i=Random(0, j-1);
	i=NV[i];
	free(NV);
	return i;
}
// generate a maze
Maze InitMaze(int Nx, int Ny)
{
	Maze M={0};
	int n;
	if ((M.down=calloc(sizeof(char),Nx*Ny))==NULL)
		return M;
	if ((M.right=calloc(sizeof(char),Nx*Ny))==NULL)
	{
		free(M.down);
		M.down=NULL;
		return M;
	}
	if ((M.visited=calloc(sizeof(char),Nx*Ny))==NULL)
	{
		free(M.visited);
		M.visited=NULL;
		return M;
	}
	M.visited[0]=1; // mark the top left corner as visited
	// all memory allocated, good to go
	M.Nx=Nx;
	M.Ny=Ny;
	
	n=SelectStart(&M);
	while (n>=0)
	{
		RandomWalk(&M, n);
		n=SelectStart(&M);
	}
	
	return M;	
}

// place the robot in the maze
curiosity InitCuriosity()
{
	curiosity C;
	int o;
	C.n=0;
	
	o=Random(0,3);
	switch(o)
	{
		case 0:
			C.o='u';
			break;
		case 1:
			C.o='r';
			break;
		case 2:
			C.o='d';
			break;
		case 3:
		default:
			C.o='l';
			break;
	}
	return C;
}

// drive curiosity around
#define NAV_SUCCESS 	0
#define NAV_CRASH 		1
#define NAV_LOST 		2
#define NAV_ILLEGAL		3
#define NAV_HITNRUN		4
route NavigateMaze(curiosity *C, char *program, Maze *M, int *success)
{
	char *p;
	char *O="urdl";
	int i, r, c, n, k=0;
	route R;
	
	R.N=strlen(program);
	R.n=malloc(R.N*sizeof(int));
	R.pd=calloc(R.N,sizeof(int));
	p=program;
	// set orientation index i
	switch(C->o)
	{
		case 'u':
			i=0;
			break;
		case 'r':
			i=1;
			break;
		case 'd':
			i=2;
			break;
		case 'l':
			i=3;
			break;
		default:
			i=0;
	}
	while (*p)
	{
		R.n[k]=C->n;
		if (C->n==M->Nx*M->Ny-1)
		{
			(*success)=NAV_HITNRUN;// curiosity drives right over the delicate life form
			goto route_return;
		}
		// rotate
		switch (*p)
		{
			case 'l': // rotate 90 degrees counter clockwise
				i--;
				switch(i)
				{
					case -1: //u-l
						R.pd[k]=R_up_left;
						break;
					case 0: //r-u
						R.pd[k]=R_right_up;
						break;
					case 1: //d-r
						R.pd[k]=R_down_right;
						break;
					case 2: //l-d
						R.pd[k]=R_left_down;
						break;
				}						
				break;
			case 'r': // rotate 90 degrees clockwise
				i++;
				switch(i)
				{
					case 1: //u-r
						R.pd[k]=R_up_right;
						break;
					case 2: //r-d
						R.pd[k]=R_right_down;
						break;
					case 3: //d-l
						R.pd[k]=R_down_left;
						break;
					case 4: //l-u
						R.pd[k]=R_left_up;
						break;
				}	
				break;
			case 'b': // rotate 180 degrees
				i+=2;
				switch(i)
				{
					case 2: //u-d
						R.pd[k]=R_up_down;
						break;
					case 3: //r-l
						R.pd[k]=R_right_left;
						break;
					case 4: //d-u
						R.pd[k]=R_down_up;
						break;
					case 5: //l-r
						R.pd[k]=R_left_right;
						break;
				}	
				break;
			case 'f': // no rotation
				switch(i)
				{
					case 0: //u-u
						R.pd[k]=R_vertical;
						break;
					case 1: //r-r
						R.pd[k]=R_horizontal;
						break;
					case 2: //d-d
						R.pd[k]=R_vertical;
						break;
					case 3: //l-l
						R.pd[k]=R_horizontal;
						break;
				}	
				break;
			default: // unknown instruction
				(*success)=NAV_ILLEGAL;
				goto route_return;
				break;
		}
		k++;
		if (i<0)
			i+=4;
		i=i%4;
		C->o=O[i];
		// move
		r=ROW(C->n,M);
		c=COL(C->n,M);
		switch (C->o)
		{
			case 'u':
				r--;
				if (r<0) // crash!
				{
					(*success)=NAV_CRASH;
					goto route_return;
				}
				n=INDEX(r,c,M);
				if (M->down[n]==0)
				{
					(*success)=1;
					goto route_return;
				} // crash!
				C->n=n;
				break;
			case 'd':
				if (r==M->Ny-1) // crash!
				{
					(*success)=NAV_CRASH;
					goto route_return;
				}
				if (M->down[C->n]==0)
				{
					(*success)=NAV_CRASH;
					goto route_return;
				} // crash!
				r++;				
				C->n=INDEX(r,c,M);
				break;
			case 'l':
				c--;
				if (c<0) // crash!
				{
					(*success)=NAV_CRASH;
					goto route_return;
				}
				n=INDEX(r,c,M);
				if (M->right[n]==0)
				{
					(*success)=NAV_CRASH;
					goto route_return;
				} // crash!
				C->n=n;
				break;
			case 'r':
				if (c==M->Nx-1) // crash!
				{
					(*success)=NAV_CRASH;
					goto route_return;
				}
				if (M->right[C->n]==0)
				{
					(*success)=NAV_CRASH;
					goto route_return;
				} // crash!
				c++;				
				C->n=INDEX(r,c,M);
				break;
			default:
				break;
		}
		p++; // advance in the program
	}
	// test target:
	if (C->n==M->Nx*M->Ny-1)
		(*success)=NAV_SUCCESS;// robot finds kitten!
	else
		(*success)=NAV_LOST;// curiosity got lost in the maze
route_return:
	R.N=k;
	return R;
} 		

// get a program
#define MAXPLEN 255
char *GetProgram()
{
	char *prog, *r, *p;
	prog=calloc((MAXPLEN+1), sizeof(char));
	printf("Curiosity Navigational Prompt:\n> ");
	if (fgets(prog,MAXPLEN,stdin)==NULL)
		prog[0]='\0';
	r=prog;
	p=prog;
	// strip endlines and numbers and whatever
	while (*p)
	{
		*r=*p;
		if(isalpha(*p))
			r++;
		p++;
	}
	*r=*p;
	
	return prog;
}

void FreeMaze(Maze *M)
{
	if (M->down)
		free(M->down);
	if (M->right)
		free(M->right);
	if (M->visited)
		free(M->visited);
	M->down=NULL;
	M->right=NULL;
	M->visited=NULL;
	M->Nx=0;
	M->Ny=0;	
}

void FreeRoute(route *R)
{
	if (R->n)
		free(R->n);
	if (R->pd)
		free(R->pd);
	R->n=NULL;
	R->pd=NULL;
	R->N=0;
}

char *WriteStr(char *p, size_t msize, const char *format, ...)
{
    va_list ap;
    size_t l=0;
    va_start(ap, format);
    
    vsnprintf(p, msize,format, ap);
    while((*p)&&(l<msize))
    {
		p++;
		l++;
	}
	return p;
} 

int RouteDirIndex(int n, route *R)
{
	int i=0;
	if (!R)
		return 14;
	while (i<R->N)
	{
		if (R->n[i]==n)
			return R->pd[i];
		i++;
	}
	return 14;
}
 
char *RouteDirString(int n, route *R) 
{
	int i;
	i=RouteDirIndex(n,R);
	return paths[i];
}

#define CELLLEN 20
char *WriteMaze(Maze *M, curiosity *C, route *R)
{
	int i=0,j=0, n;
	int B;
	int r,c;
	char *cc;
	char *out, *p;
	int Len;
	
	Len=(M->Nx+3)*(M->Ny*2+2)*CELLLEN;
	
	if((out=calloc(Len,sizeof(char)))==NULL)
		return NULL;
	
	r=ROW(C->n,M);
	c=COL(C->n,M);
	switch (C->o)
	{
		case 'u':
			cc=Cu;
			break;
		case 'd':
			cc=Cd;
			break;
		case 'r':
			cc=Cr;
			break;
		case 'l':
		default:
			cc=Cl;
			break;
	}
	p=out;		
	// top fat border
	p=WriteStr(p,CELLLEN, " ");
	for (i=0;i<M->Nx;i++)
		p=WriteStr(p,CELLLEN, " %c  ", 'A'+i);
	p=WriteStr(p,CELLLEN, "\n");
	p=WriteStr(p,CELLLEN, outerborders[0]);
	for (i=0;i<M->Nx-1;i++)
	{
		if (M->right[INDEX(j,i,M)])
			p=WriteStr(p,CELLLEN, outerborders[1]);
		else
			p=WriteStr(p,CELLLEN, outerborders[2]);
	}
	p=WriteStr(p,CELLLEN,outerborders[3]);
	
	// body			
	for (j=0;j<M->Ny-1;j++)
	{
		// left fat border
		p=WriteStr(p,CELLLEN,"\n");
		p=WriteStr(p,CELLLEN,outerborders[4]);
		for (i=0;i<M->Nx-1;i++)
		{
			n=INDEX(j,i,M);
			if ((i==c)&&(j==r))
				p=WriteStr(p,CELLLEN," %s ",cc);
			else
				p=WriteStr(p,CELLLEN," %s ",RouteDirString(n, R));			
			
			if (M->right[INDEX(j,i,M)])
				p=WriteStr(p,CELLLEN," ");
			else
				p=WriteStr(p,CELLLEN,innerborders[16]);
		}
		n=INDEX(j,i,M);
		// right fat border
		p=WriteStr(p,CELLLEN," %s ",RouteDirString(n, R));
		p=WriteStr(p,CELLLEN,outerborders[4]);
		p=WriteStr(p,CELLLEN, " %d  ", j+1);
		
		// left fat border
		p=WriteStr(p,CELLLEN,"\n");
		if (M->down[INDEX(j,0,M)])
			p=WriteStr(p,CELLLEN,outerborders[4]);
		else
			p=WriteStr(p,CELLLEN,outerborders[5]);
				
		for (i=0;i<M->Nx-1;i++)
		{
			B=(M->down[INDEX(j,i,M)]==1);
			B|=2*(M->right[INDEX(j,i,M)]==1);
			B|=4*(M->down[INDEX(j,i+1,M)]==1);
			B|=8*(M->right[INDEX(j+1,i,M)]==1);
			p=WriteStr(p,CELLLEN,innerborders[B]);
		}
		// right fat border
		if (M->down[INDEX(j,i,M)]==1)
		{
			p=WriteStr(p,CELLLEN,"   ");
			p=WriteStr(p,CELLLEN,outerborders[4]);
		}
		else
			p=WriteStr(p,CELLLEN,outerborders[6]);
				
	}
	// left fat border
	p=WriteStr(p,CELLLEN,"\n");
	p=WriteStr(p,CELLLEN,outerborders[4]);
	
	for (i=0;i<M->Nx-1;i++)
	{
		n=INDEX(j,i,M);
		if ((i==c)&&(j==r))
			p=WriteStr(p,CELLLEN," %s ",cc);
		else
			p=WriteStr(p,CELLLEN," %s ",RouteDirString(n, R));		
		if (M->right[INDEX(j,i,M)])
			p=WriteStr(p,CELLLEN," ");
		else
			p=WriteStr(p,CELLLEN,innerborders[16]);
	}
	if (C->n==M->Nx*M->Ny-1)
		p=WriteStr(p,CELLLEN," X ");
	else
		p=WriteStr(p,CELLLEN," O ");
	p=WriteStr(p,CELLLEN,outerborders[4]);
	p=WriteStr(p,CELLLEN, " %d  ", j+1);
	p=WriteStr(p,CELLLEN,"\n");
	p=WriteStr(p,CELLLEN,outerborders[7]);
	
	for (i=0;i<M->Nx-1;i++)
	{
		if (M->right[INDEX(j,i,M)])
			p=WriteStr(p,CELLLEN,outerborders[1]);
		else
			p=WriteStr(p,CELLLEN,outerborders[8]);
	}
	p=WriteStr(p,CELLLEN,outerborders[9]);
	p=WriteStr(p,CELLLEN,"\n");
	return out;
}		

void PrintMaze(Maze *M, curiosity *C, route *R)
{
	char *out;
	out=WriteMaze(M,C,R);
	printf("%s\n",out);
	free(out);
}

char *FindPause(char *str, int *P)
{
	char *p, *q;
	(*P)=0;
	p=str;
	while (*p)
	{
		if (strncmp(p,"<P",2)==0)
		{
			q=p+2;
			while ((*q>='0')&&(*q<='9'))
				q++;
			if (*q=='>')
			{
				*p='\0';
				(*P)=atoi(p+2);
				q++;
				*p='\0';
				return q;
			}
		}
		p++;
	}
	return NULL;
}
	

void PrintMessage(char *message)
{
	char *chunked;
	char *p, *q;
	int pause=0;
	chunked=malloc((strlen(message)+1)*sizeof(char));
	p=chunked;
	q=message;
	/* first just copy */
	while (*q)
	{
		*p=*q;
		p++;
		q++;
	}
	*p=*q;
	p=chunked;
	q=p;
	while (q)
	{
		q=FindPause(p, &pause);
		printf("%s",p);
		fflush(stdout);
		sleep(pause);
		p=q;
	}
	free(chunked);
}

void TheEnd(int r)
{
	switch (r)
	{
		case NAV_SUCCESS:
			printf("Curiosity Acquired Target!\n");
			PrintMessage(Msuccess[Random(0,NMsuccess-1)]);
			sleep(2);
			printf("NASA thanks you for your services!\n");
			break;
		case NAV_CRASH:
			PrintMessage(Mcrash[Random(0,NMcrash-1)]);
			break;
		case NAV_LOST:
			PrintMessage(Mlost[Random(0,NMlost-1)]);
			break;
		case NAV_ILLEGAL:
			PrintMessage(Millegal[Random(0,NMillegal-1)]);
			break;
		default:
			PrintMessage(Mhitnrun[Random(0,NMhitnrun-1)]);
			break;
	}
}

void Curiosity(int Nx, int Ny)
{
	Maze M;
	curiosity C;
	route R;
	char *prog;
	int r;
	char *cc;
	InitRandom();
	M=InitMaze(Nx, Ny);
	C=InitCuriosity();
	switch (C.o)
	{
		case 'u':
			cc=Cu;
			break;
		case 'd':
			cc=Cd;
			break;
		case 'r':
			cc=Cr;
			break;
		//case 'l':
		default:
			cc=Cl;
			break;
	}
	printf("%s",Intro);
	printf("-- Press Enter To Continue --");
	getchar();
	printf(Situation,cc);
	PrintMaze(&M,&C,NULL);
	printf("-- Press Enter To Continue --");
	getchar();
	printf(Guide,cc);
	printf("-- Press Enter To Continue --");
	getchar();
	PrintMaze(&M,&C,NULL);
	
	prog=GetProgram();
	printf("Please stand by as the program is\nuploaded...\n");
	sleep(1);
	printf("Curiosity acknowledges receipt.\ncompiling...\n");	
	sleep(1);
	printf("Executing, please stand by...\n");	
	R=NavigateMaze(&C, prog, &M, &r);
	free(prog);
	sleep(2);
	PrintMaze(&M,&C,&R);
	TheEnd(r);
	FreeMaze(&M);
	FreeRoute(&R);
}

