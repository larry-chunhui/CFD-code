
//��ʽ������Alt+ F8
//����߽綼�޷��������ɣ����д������ɵ�������
//trinum[i]���������trinum[],�ڵڶ���ѭ�����Ѿ���ȷ�ˣ�
//nodenum �ڵ����
//ntri����������
struct Node //Ϊʲô�ҵ���������������Ӻ��� ��������ͬ�أ�
{   float x;
    float y;
}   node[10000],mode[1000][1000];//�ӵ�һ��node1��ʼ����,mode ǰ���Ǹ��Ǹ�λ

struct lnode //�ڵ�����
{   float x,y;
    int num;

    struct lnode *next;
};

struct Ltri
	{   int n1,n2,n3;  //���������εĽڵ�
    int t1,t2,t3;
}   ltri[4];

struct tri//����������
{int n1,n2,n3;  //���������εĽڵ�
 int t1,t2,t3;  //���������εı�ţ�������ʱ��˳��
 int num;//�����α��
 struct tri *next;
	 
};

#define NULL 0
#define LEN sizeof (struct lnode)
#define len sizeof (struct tri)
#include<iostream>
#include<stdio.h>
#include<stdlib.h> 
#include <math.h>
const char file_name[50] = "d:\\����ͨ��.txt";  

int main(int argc, char *argv[])
{ struct lnode *head; 
  float d(float x1,float y1,float x2,float y2);

  void readfile( Node(&node)[100000],int &a, FILE *fp1);     //���ļ�������򣬲����ɳ�ʼ����
  struct lnode *nodeshengcheng(Node(&node)[100000],int &m);//�����ڵ����� 
  struct lnode *del(struct lnode *head, int i);  //ɾ���ڵ������б��Ϊi�Ľڵ�
  struct lnode *delrongru(struct lnode *head);   //ɾ���ڵ��������ظ��ĵ�
  void print(struct lnode *head);                //��ʾ���ڵ��������ڼ��  
  Node node[100000]; 
  
  void dengfenbianjie(Node(&node)[100000],Node(&mode)[1000][1000],int &m,int &h,float &b,int NN);

  struct tri *thead;                             //��ʼ������ͷָ��
  struct tri *tcreat (Ltri(&ltri)[4],int &m);     //��������������

  struct tri *deltri(struct tri *thead, int i);  //ɾ�������������б��Ϊi��������
  struct tri *xuanzhuan(struct tri *thead, int i);  //��ת�����������б��Ϊi��������
  struct tri *deltriall(struct tri *thead,struct tri *p );//ֻɾ��ĳһ��������
  struct tri *deltribianjie(struct tri *thead, int &m);  //ɾ���߿��ĸ������������
  struct tri *deltrisuoyou( Node(&node)[100000],struct lnode *head,struct tri *thead,int &m);  //ɾ���������е�������
  struct tri *bianjietri(struct tri *thead,Node(&node)[100000],int &m,int &h);//���ñ߽�ڵ�����������
  void tprint(struct tri *thead);                 //��ʾ���������������ڼ��
  int  count (struct tri *thead);                  //ͳ������������Ŀ�����ڶ������ӵ������α��
  int inside(Node(&node)[100000],struct tri *thead,struct lnode *head,int (&trinum)[1000],int &ntri,int &nodenum) ;//�������κͽڵ�
  void newd(struct tri *thead,int (&trinum)[1000],int &ntri,int &nodenum);//delauny���������������ڽڵ�
  void newi(struct tri *thead,int (&trinum)[1000],int &ntri,int &nodenum);//delauny����һ���������ڽڵ�
  void intri(Node(&node)[100000],struct tri *thead,int trinum[],int nodenum,int n);                            //�����������������Բ������a���ҳ����������β���ڵ�
  void delaunay(Node(&node)[100000],struct tri *thead);//�����������ν���Delaunay ���߲��� 
  void writefile( Node(&node)[100000],int &h,struct tri *thead, FILE *fp2);//д�ļ����ڵ�����ʦ�����
  int online(Node(&node)[100000],int a,int b,int c);//�ж�c�Ƿ���ֱ��ab��,����ֵΪ1��ʾ�������ϣ�0��ʾ�㲻������
  int satisfy(Node(&node)[100000],struct tri *p,float &b);//�ж�������p�Ƿ����������߶��ӽ��ڱ߳�b��Ҫ�����㷵��1�������㷵��0
  void insert(struct tri *thead,Node(&node)[100000],float b,int &h);//��ʽ���ʼ�����ڲ�����㣬���ɵȱ�������ֱ�����������ζ��㹻С
  void sinsert (struct tri *thead,Node(&node)[100000],float &b,int &h);//2.������ܴ���������ڲ���㣬����ʣ��������
   void deinsert(struct tri *thead,Node(&node)[100000],float &b,int &h);//�����������ģ��
  //int  ornot(Node(&node)[100000],struct tri *p,float x1,float y1);//�жϵ�(x,y)�Ƿ���tri p �ڣ�����=1������=0
  float s(float x[3],float y[3]);//��������p�����
  void youhua(struct tri *thead,Node(&node)[100000],float &b);//ɾ��̫С�������Σ��������ڵ�ϲ���һ������������������
  void lapace(struct tri *thead,Node(&node)[100000],int &hh,int &h);//lapace �Ż�����
  float  mintri(struct tri *thead,Node(&node)[100000],int i);//��������i��̱߳�
    float  maxtri(struct tri *thead,Node(&node)[100000],int i);//��������i��߳�
  //ר�����ڵ��Գ���ĺ���
    void WRitefile( Node(&node)[100000],int &h,struct tri *thead, FILE *fp3);//д�ļ����ڵ�����ʦ�����
	  void qualitytri(Node(&node)[100000],struct tri *thead);
  int NN=3000;
  int i=1;//i�������ڼ���
  int m=1,h=0,hh;//m�洢��ʼ�ļ����ж��ٽڵ㣬h�洢���ж��ٸ��ڵ�,hh�洢���ж��ٸ��߽�ڵ�
  float b;//b��ʾ�ȱ������α߳�
  int trinum[1000];
  int ntri=0;
  int nodenum=0;
  int f; //�ռ�inside�ķ���ֵ
  FILE *fp1,*fp2,*fp3;
 
	  readfile(node,m,fp1);    //�����ļ�
	  head =nodeshengcheng(node,m);                           //�����߽�ڵ�����
	  head=delrongru(head);//���벻������ˣ���Ϊ
	  thead=tcreat (ltri,m);//������ʼ���������ε�����
	  
QQ: if (head!=NULL)	  
	{ 
		
		f=inside(node,thead,head,trinum,ntri,nodenum);
		
		if(f==1)
		{      newd(thead,trinum,ntri,nodenum);
		       delaunay(node,thead);
		
		
		}
		else if(f==0)
		{      newi(thead,trinum,ntri,nodenum);
	           delaunay(node,thead);
		}
		
		head=del(head,nodenum);
		
		goto QQ;
		
	}

	thead=deltribianjie(thead, m);//ɾ����߽��ϵ�������
	delaunay(node,thead);
    thead=deltrisuoyou( node,head,thead,m); //ɾ����ͨ�������������



	dengfenbianjie(node,mode,m,h,b,NN);
	thead=bianjietri(thead,node,m,h);//���ñ߽�ڵ�����������
	delaunay(node,thead);//�������������
	for(i=m+1;i<=m+h;i++)
	{	       
		if (online(node,1,2,i)==1)
		{	node[i].x=0;
	    	node[i].y=0;
		}
	}
	
	h=h+m;

	hh=h;
	insert(thead,node,b,h);

	sinsert (thead,node,b,h);//2.������ܴ���������ڲ���㣬����ʣ��������
	printf("     woqunimei\n");
			writefile( node,h,thead, fp2);		
	//���������Ż�������400��������ʱ�޷����У�����
//	youhua(thead,node,b);
		printf("woqunimei\n");
	lapace(thead,node,hh,h);
 

	
			
			qualitytri(node,thead);
			
			
			
			



return 0;
}

float d(float x1,float y1,float x2,float y2)//���������뺯��
{
	
	float da;
	da=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	return da;
 }


void readfile(Node(&node)[100000],int &m,FILE *fp1)//������򣬲����ɳ�ʼ����
{   FILE *fp;
  	int i=0,j=0; 
	int cout=0;//�����ɱ߽�ڵ�
	float a[1000]; //�����б߽�ڵ����ݵ������һλ������������ת
	float temp; 
	float xmax, ymax,xmin,ymin;
	      fp1=fopen("d:\\��ʼ�������нڵ�.txt","w");     //�������������ʼ�������нڵ�.txt
	if((fp=fopen(file_name, "rb")) == NULL)
	{     printf("��ȷ���ļ�(%s)�Ƿ����!\n", file_name);    
	      exit(1);    
	}   
	
	while(1==fscanf(fp, "%f", &temp)) 
	{   a[i]=temp ;
	    i++;
	} 
	    m=(i-1)/2;
	for(i=1;i<=m;i++)
	{   node[i].x=a[2*i-1];
	    node[i].y=a[2*i];
	}




    xmax=xmin=node[1].x;
	ymax=ymin=node[1].y;
	for (i=2;i<=m;i++)
	{   if (xmax<node[i].x)
	        xmax=node[i].x;
	    if (xmin>node[i].x)
		    xmin=node[i].x;
	    if (ymax<node[i].y)
		   ymax=node[i].y;
	    if (ymin>node[i].y)
		   ymin=node[i].y;
	}
	
	xmax=xmax+35;
	ymax=ymax+35;
	xmin=xmin-35;
	ymin=ymin-35;
	node[m+1].x=xmin;
	node[m+1].y=ymax;
	node[m+2].x=xmin;
	node[m+2].y=ymin;
	node[m+3].x=xmax;
	node[m+3].y=ymin;
	node[m+4].x=xmax;
	node[m+4].y=ymax;
	node[m+5].x=xmin;
	node[m+5].y=ymax;
	fprintf(fp1,"%d\n",m+5);
	for(i=1;i<=m+5;i++)
		fprintf(fp1,"%f  %f\n",node[i].x,node[i].y);
	fclose (fp);


}

struct lnode *nodeshengcheng (Node(&node)[100000],int &m)//�����߽�ڵ�����
{	int i;
	struct lnode *head;//����˫�������е�ָ��ڵ�
	struct lnode *p1;
	struct lnode *p2;
	p1=p2=(struct lnode *)malloc(LEN);// �����������Լ�������ȥ�Ĵ���
	for (i=1;i<=m;i++)
	{
		p1->x=node[i].x;
	    p1->y=node[i].y;
	    p1->num=i;
	 
	if (i==1)
	{   head =p1;
	}
	else p2->next =p1;
	     p2=p1;
	     p1=(struct lnode *)malloc(LEN);
	}
	     p2->next =NULL;
		 return(head);

}

void print(struct lnode *head)//���������
{        struct lnode *p;
	     printf("lnode links are \n ");  //�������
     	 p=head ;
	     if (head !=NULL)
		 do{
	     	printf("i=%d  %f   %f\n",p->num,p->x,p->y);
			p=p->next;
			
		}
		while (p!=NULL);
}



struct lnode *del(struct lnode *head, int i)
{struct lnode *p1,*p2;
    if (head==NULL) {printf("\nlist null!\n");return head;}
	 p1=head;

	while (i!=p1->num && p1->next !=NULL)
	{  p2=p1;p1=p1->next;}
	if (i==p1->num)
	{if (p1==head)  
	{head=p1->next;
	 return (head);
	}
	else {p2->next=p1->next;return (head);}
	 printf("delete:%d\n",i);
	}

	else printf("%d not found !\n",i);
	return (head);
}

 struct lnode *delrongru(struct lnode *head)//ɾ���ڵ��������ظ��ĵ�
 { 
	 struct lnode *p,*q;
	 for(p=head;p!=NULL;p=p->next)
	 { 
		 for(q=p->next;q!=NULL;q=q->next)
		 { 
			 if (p->x==q->x && p->y==q->y)
		      {
			   head=del(head,q->num);
		       
			  }
		 
		 }
	 }
	 return (head);

}

struct tri *tcreat (Ltri(&ltri)[4],int &m)                                 //��������������
{   struct tri *p11,*p22,*thead;
    int i=1;
    ltri[1].n1=m+1;ltri[1].n2=m+2;ltri[1].n3=m+3;
	ltri[1].t1=0;ltri[1].t2=0;ltri[1].t3=2;
	ltri[2].n1=m+3;ltri[2].n2=m+4;ltri[2].n3=m+1;
	ltri[2].t1=0;ltri[2].t2=0;ltri[2].t3=1;

		p11=p22=(struct tri *)malloc(len);// �����������Լ�������ȥ�Ĵ���
	for (i=1;i<=2;i++)
	{
		p11->n1=ltri[i].n1;
		p11->n2=ltri[i].n2;
		p11->n3=ltri[i].n3;
		p11->t1=ltri[i].t1;
		p11->t2=ltri[i].t2;
		p11->t3=ltri[i].t3;
		p11->num=i;
      
	 
	if (i==1)
	{   thead =p11;
	}
	else p22->next =p11;
	     p22=p11;
	     p11=(struct tri *)malloc(len);
	}
	     p22->next =NULL;
		 return(thead);
}
float s(float x[3],float y[3])//��������p�����
{      
	
	   float a,b,c,e,f;
	   
	   a=d(x[1],y[1],x[2],y[2]);
	   b=d(x[3],y[3],x[2],y[2]);
	   c=d(x[1],y[1],x[3],y[3]);
	   e=0.5*(a+b+c);
	   f=sqrt(e*(e-a)*(e-b)*(e-c));
	   return f;
	   
}

struct tri *xuanzhuan(struct tri *thead, int i) //��ת�����������б��Ϊi��������
{   struct tri *p;
int k;
for(p=thead;p->num!=i;p=p->next)
{};

k=p->n1;       
p->n1=p->n2;
p->n2=p->n3;
p->n3=k;
k=p->t1;       
p->t1=p->t2;
p->t2=p->t3;
p->t3=k;
return (thead);
 }


void tprint(struct tri *thead)
{        struct tri *p;
         printf("tri links are\n ");  //�������
     	 p=thead ;
	     if (thead !=NULL)
		 do{
			printf("trinum=%d\n",p->num);
	     	printf("%d  %d  %d\n",p->n1,p->n2,p->n3);
	        printf("%d  %d  %d\n",p->t1,p->t2,p->t3);
			p=p->next;
			}
		while (p!=NULL);
}


int max(struct tri *thead)//��������������������
{   int i=0;
struct tri *p;
for(p=thead;p!=NULL;p=p->next)
{
	if (i<p->num)
	{i=p->num;}
}
return i;
}



 int count (struct tri *thead)
 {int i=0;
  struct tri *p;
  for(p=thead;p!=NULL;i++,p=p->next)
  {}
  return i;
 
 }

 int inside(Node(&node)[100000],struct tri *thead,struct lnode *head,int (&trinum)[1000],int &ntri,int &nodenum)  //�жϵ������������������Բ�ڣ������������������������Բ�У�����������������ţ�
	                                                             //��ֻ���������������Բ�ڣ��򷵻�ֵΪ1����ֻ��һ���������У��򷵻�ֵΪ0
{          struct tri *p;                               
           struct lnode *q;                                      //numnode �ڵ���ţ�node[numnode].x,node[numnode].y                         
		          q=head;
           float  x1,x2,x3,y1,y2,y3;  
           float  x[1000] ;  //Բ�ĺ�����
           float  y[1000];  //Բ��������
		   float  r[1000];      //Բ�İ뾶
		   int    i=0,m=0,j=1;   //i���ڼ�����m��¼���ж��ٸ�������
		   int    n;            //n��¼�ж������������Բ�����ڵ�
           
           void intri(Node(&node)[100000],struct tri *p,int trinum[],int nodenum,int n);//����һ���Ǹ������������������Բ�����ڵ�ĺ���
           for (i=1,p=thead;p!=NULL;p=p->next,i++)//������������Բ�ģ��뾶
		      { x1=node[p->n1].x; 
                x2=node[p->n2].x;   
                x3=node[p->n3].x; 
                y1=node[p->n1].y; 
                y2=node[p->n2].y; 
                y3=node[p->n3].y;   
				x[i]=((y2-y1)*(y3*y3-y1*y1+x3*x3-x1*x1)-(y3-y1)*(y2*y2-y1*y1+x2*x2-x1*x1))/(2*(x3-x1)*(y2-y1)-2*((x2-x1)*(y3-y1)));  
                y[i]=((x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1)-(x3-x1)*(x2*x2-x1*x1+y2*y2-y1*y1))/(2*(y3-y1)*(x2-x1)-2*((y2-y1)*(x3-x1)));  
                r[i]=sqrt((x1-x[i])*(x1-x[i])+(y1-y[i])*(y1-y[i]));

		   }
		   ntri=i-1;                                                            //i=3,m=2,m������������������������
		  	   
				

p=thead;	//��pָ�������������׵�ַ	   
LOOP:for (i=1;i<=ntri;i++)              //����m�������Σ����������������ҳ����нڵ��
{
if ((q->x-x[i])*(q->x-x[i])+(q->y-y[i])*(q->y-y[i])-r[i]*r[i]<0)//qָ���׸��ڵ� ���������Բ��
				{      trinum[j]=i;            //������Բ�����ýڵ�����������trinum[i]������i
				       j++;
				}
}



	             n=(j-1);                            //��¼�ڽ�Բ�����ýڵ������������
				 if (n>2)                           //����������������Σ��򷵻�1
				      {nodenum=q->num;
				       intri(node,thead,trinum,nodenum,n);

				       return 1;}
				 if (n==2)                           //����������������Σ��򷵻�1
				      {nodenum=q->num;
				      return 1;}
                      else if(n==1)
					         {
					          if (q->next==NULL)
							       {nodenum=q->num;  //�����ֻ��һ���������ϣ��򷵻�0   
							        return 0;}  
					          else 
							  {   q=q->next; 
								  goto LOOP;//�����ж���һ���ڵ��Ƿ��������������������Բ��
							  }
							  
								  
					  }
} 


void newd(struct tri *thead,int (&trinum)[1000],int &ntri,int &nodenum)
 { struct tri *d1,*d2,*p,*p0;
   int n,m;
   ntri= max (thead);
   for(p=thead;p->num!=trinum[1];p=p->next)//��d1ָ���һ��������
   {};
          d1=p;

   for(p=thead;p->num!=trinum[2];p=p->next)//��d2ָ��ڶ���������
   {};
          d2=p;
		  do                       //��d1��������ת���������λ��
		  {n=d1->n1;       
		  d1->n1=d1->n2;
		  d1->n2=d1->n3;
		  d1->n3=n;
		  m=d1->t1;       
		  d1->t1=d1->t2;
		  d1->t2=d1->t3;
		  d1->t3=m;
		  }while(d1->n2==d2->n1 || d1->n2==d2->n2  ||  d1->n2==d2->n3);

		  do                       //��d2��������ת���������λ��
		  {n=d2->n1;       
		  d2->n1=d2->n2;
		  d2->n2=d2->n3;
		  d2->n3=n;
		  m=d2->t1;       
		  d2->t1=d2->t2;
		  d2->t2=d2->t3;
		  d2->t3=m;
		  }while(d2->n2==d1->n1 || d2->n2==d1->n2  ||  d2->n2==d1->n3);
 
for(p=thead ;p->next!=NULL;p=p->next)//��ָ��pָ���β֮��Ŀ�����
   {};


    
   {      p0=(struct tri *)malloc (len);
	      p0->num=ntri+1;          //�ڶ�β�����һ��������
          p0->n1=nodenum;
		  p0->n2=d1->n2;
		  p0->n3=d2->n1;
		  p0->t1=d1->num;
		  p0->t2=d1->t2;
		  p0->t3=d2->num;
		  p->next=p0;
		  p0->next=NULL;
   } 
    


  

          p0=(struct tri *)malloc (len);
          p=p->next;
		 {p0->num=ntri+2;          //�ڶ�β����ڶ���������
		  p0->n1=nodenum;
		  p0->n2=d2->n2;
		  p0->n3=d1->n1;
		  p0->t1=d2->num;
		  p0->t2=d2->t2;
		  p0->t3=d1->num;
		 }p->next=p0;
		  p0->next=NULL;



         {d1->n3=nodenum;         //�滻ԭ����һ��������
		  d1->t2=ntri+1;
		  d1->t3=ntri+2;
		 }

		 {d2->n3=nodenum;        //�滻ԭ���ڶ���������
		  d2->t2=ntri+2;
		  d2->t3=ntri+1;
		 }
		 


		    for(p=thead;p->num!=ntri+1;p=p->next)//��d1ָ�����Ϊntri+1��������
			{};
            d1=p;

            for(p=thead;p->num!=ntri+2;p=p->next)//��d2ָ�����Ϊntri+2��������
			{};   
            d2=p;


		n=d1->t2;//�ҵ���һ�������εڶ������ڱ�t2����Ӧ����������ţ����ı���t2��ֵ������Ϊ��ʹԭ��������������Χ���ĸ�������t1,t2,t3��ȷ��
		if (n!=0)
		{
			for(p=thead;p->num!=n;p=p->next)
			{};
	         if (p->t1==trinum[1])          p->t1=(ntri+1);//���ˣ����������1ĳһ��ָ������������
		     else if (p->t2==trinum[1])     p->t2=(ntri+1);
		     else if (p->t3==trinum[1])     p->t3=(ntri+1);
		}

		   
		                 	   

		n=d2->t2;
		if (n!=0)
		{   
			for(p=thead;p->num!=n;p=p->next)
			{};
		    if (p->t1==trinum[2])           p->t1=(ntri+2);
		    else if (p->t2==trinum[2])      p->t2=(ntri+2);
		    else if (p->t3==trinum[2])      p->t3=(ntri+2); 
		}


		
 }


void newi(struct tri *thead,int (&trinum)[1000],int &ntri,int &nodenum)//���������в���һ����
{  
	struct tri *d,*p,*p0,*p1;
    int n;

   for(p=thead;p->num!=trinum[1];p=p->next)//��dָ�������������
   {};
          d=p;


		  for(p=thead ;p->next!=NULL;p=p->next)//��ָ��pָ���β������
   {};    p0=(struct tri *)malloc (len);
   {      p0->num=ntri+1;          //�ڶ�β�����һ��������
          p0->n1=nodenum;
		  p0->n2=d->n2;
		  p0->n3=d->n3;
		  p0->t1=d->num;
		  p0->t2=d->t2;
		  p0->t3=ntri+2;
		  p->next=p0;
		  p0->next=NULL;

   }      p1=(struct tri *)malloc (len);
          p=p->next;
		 {p1->num=ntri+2;          //�ڶ�β����ڶ���������
		  p1->n1=d->n3;
		  p1->n2=d->n1;
		  p1->n3=nodenum;
		  p1->t1=d->t3;
		  p1->t2=d->num;
		  p1->t3=ntri+1;
		  p->next=p1;
		  p1->next=NULL;
		 }
		 
		 {d->n3=nodenum;         //�滻ԭ��������
		  d->t2=ntri+1;
		  d->t3=ntri+2;
		 }

		n=p0->t2;//�ҵ���һ�������εڶ������ڱ�t2����Ӧ����������ţ����ı���t2��ֵ������Ϊ��ʹԭ��������������Χ���ĸ�������t1,t2,t3��ȷ��
		if (n!=0)
		{
			for(p=thead;p->num!=n;p=p->next)
			{};
	         if (p->t1==trinum[1])          p->t1=(ntri+1);
		     else if (p->t2==trinum[1])     p->t2=(ntri+1);
		     else if (p->t3==trinum[1])     p->t3=(ntri+1);
		}

		n=p1->t1;//�ҵ���һ�������εڶ������ڱ�t2����Ӧ����������ţ����ı���t2��ֵ������Ϊ��ʹԭ��������������Χ���ĸ�������t1,t2,t3��ȷ��
		if (n!=0)
		{
			for(p=thead;p->num!=n;p=p->next)
			{};
	         if (p->t1==trinum[1])          p->t1=(ntri+2);
		     else if (p->t2==trinum[1])     p->t2=(ntri+2);
		     else if (p->t3==trinum[1])     p->t3=(ntri+2);
		}
}

void intri(Node(&node)[100000],struct tri *thead,int trinum[],int nodenum,int n)//�����������������Բ������a���ҳ����������β���ڵ�
{  
   struct tri *p;
   float x[4],y[4],f[3];
   float x1;
   float y1;
   int i,j,k,a,b,a1,a2,a3;
   x1=node[nodenum].x;
   y1=node[nodenum].y;

   for(k=1;k<=n;k++)                                 //�ҵ��������Բ����nodenum��������
  {   

       
	   for (p=thead;p->num!=trinum[k];p=p->next)        //pָ��trinum[k]���������
	   { };
         a1=p->n1;a2=p->n2;a3=p->n3;//����ֱ�������㷨
         x[1]=node[a1].x;
		 y[1]=node[a1].y;
		 x[2]=node[a2].x;
		 y[2]=node[a2].y;
		 x[3]=node[a3].x;
		 y[3]=node[a3].y;
		 x[4]=node[a1].x;
		 y[4]=node[a1].y;		 

		 j=0;                                        //����ǳ���Ҫ��ÿ�������ζ�Ҫ�Ƚ�j��0
		 
		 for(i=1;i<=3;i++)                           //�ҵ�������a������nodenum
		 {   
			 if (x[i]!=x[i+1])

			 {   printf("i=%d\n",i);
				 f[i]=y[i]+(x1-x[i])/(x[i+1]-x[i])*(y[i+1]-y[i]);
			       if((y1-f[i])*(x[i+1]-x[i])>=0)
				   {
				   j++;
				   }
			 }
			 else if((x1-x[i+1])*(y[i+1]-y[i])<0)
			 {j++;}
		 }
		 
		 if (j==3)
		 { a=trinum[k];
		   break;
		 }
   }


   for(p=thead;p->num!=a;p=p->next)                    //�ҵ�n������������a���ٵ���������Ϊtrinum[2}
   {};                                                 //printf("xia yi ge xiang lin  p->num=%d\n",p->num);
   for(k=1;k<=n;k++)
   {
	   if (p->t1==trinum[k])  
	   {b=trinum[k];
	   }
       else if (p->t2==trinum[k])
	   {b=trinum[k];
	   }
       else if (p->t3==trinum[k])
	   {b=trinum[k];
	   }
   }

    trinum[1]=a;
	trinum[2]=b;

}
					 

struct tri *deltri(struct tri *thead, int i)  //ɾ�������������б��Ϊi��������					  

{struct tri *p1,*p2;
    if (thead==NULL) {printf("\nlist null!\n");return thead;}
	 p1=thead;

	while (i!=p1->num && p1->next !=NULL)
	{  p2=p1;p1=p1->next;}
	if (i==p1->num)
	{if (p1==thead)  
	{thead=p1->next;
	 return (thead);
	}
	else {p2->next=p1->next;return (thead);}
	 printf("delete:%d\n",i);
	}

	else printf("%d not found !\n",i);
	return (thead);
}

struct tri *deltriall(struct tri *thead,struct tri *p )//ֻɾ��ĳһ��������
{
	         struct tri *q;
			 int a;
	         a=p->t1;//ָ������������
			 if (a!=0)//���������⣬�������ڸ㶨��������
			 {
				 for(q=thead;q->num!=a;q=q->next)
				 {};
				 if (q->t1==p->num)      q->t1=0;
				 else if(q->t2==p->num)  q->t2=0;
				 else if(q->t3==p->num)  q->t3=0;
				 
			 }
			 
			 a=p->t2;
			 if (a!=0)
			 {
				 for(q=thead;q->num!=a;q=q->next)
				 {};
		
				 if (q->t1==p->num)      q->t1=0;
				 else if(q->t2==p->num)  q->t2=0;
				 else if(q->t3==p->num)  q->t3=0;
				
				 
			 }
			 a=p->t3;
			 if (a!=0)
			 {
				 for(q=thead;q->num!=a;q=q->next)
				 {};
				 if (q->t1==p->num)      q->t1=0;
				 else if(q->t2==p->num)  q->t2=0;
				 else if(q->t3==p->num)  q->t3=0;
				 
			 }
			 thead=deltri(thead, p->num);
			 return thead;
}

				
 struct tri *deltribianjie(struct tri *thead, int &m)  //ɾ���߿��ĸ������������,�������
 {
	 struct tri *p;
	 for (p=thead;p!=NULL;p=p->next)
	 {    
		 if (p->n1==m+1||p->n1==m+2||p->n1==m+3||p->n1==m+4 ||
			 p->n2==m+1||p->n2==m+2||p->n2==m+3||p->n2==m+4 ||
			 p->n3==m+1||p->n3==m+2||p->n3==m+3||p->n3==m+4)
		 {
			 thead=deltriall(thead,p);//������еı仯
		 }
	  
	 }

	 return thead;
 }


 
 struct tri *deltrisuoyou( Node(&node)[100000],struct lnode *head,struct tri *thead,int &m)  //ɾ���������е�������
 {       struct lnode *p,*p1,*m1,*m2;
         struct tri *t;
		 float x[3],y[3],f;
		 float i,j;
		 int a1,a2,a3,b1,b2;

		 head =nodeshengcheng(node,m); //�ָ��ڵ�����

		 for(p=head->next;node[p->num].x!=node[head->num].x ||
			 node[p->num].y!=node[head->num].y;p=p->next)
		 {};

		 p->num=head->num;//pָ����߿��һ���ڵ�,���Ϊ1
		 for(p1=p->next;node[p1->num].x!=node[head->num].x ||
			 node[p1->num].y!=node[head->num].y;p1=p1->next)
		 {};

		 p1->num=head->num;//���һ���ڵ����Ҳ��1//��ʼɾ��������
		 for(m1=p;m1->next!=NULL;m1=m1->next)//�����ڵ�,��Χ�ڵ㣡����//���������������ƻ��Ǹ㲻��
		 {   m2=m1->next;
		 x[1]=node[m1->num].x;
		 y[1]=node[m1->num].y;
		 x[2]=node[m2->num].x;
		 y[2]=node[m2->num].y;
		 for (t=thead;t!=NULL;t=t->next)//����������
		 {   
			 if( (t->n1==m1->num||t->n1==m2->num) &&(t->n2==m1->num||t->n2==m2->num))
			 {   i=node[t->n3].x;
			     j=node[t->n3].y;
			 
			          if (x[1]!=x[2])//�����������
					  {   
						  f=y[2]+(i-x[2])/(x[2]-x[1])*(y[2]-y[1]);
						  if((j-f)*(x[2]-x[1])<0)
						  {
							  thead=deltriall(thead,t);//ɾȥ���������
						  }
					  }	 
					  
					  else if((i-x[2])*(y[2]-y[1])>0)
					  {
						  thead=deltriall(thead,t);//ɾȥ���������
					  }
			 }
			 
			 else if( (t->n1==m1->num||t->n1==m2->num) &&(t->n3==m1->num||t->n3==m2->num))
			 {   i=node[t->n2].x;
			     j=node[t->n2].y;
			 if (x[1]!=x[2])//�����������
			 {   
				 f=y[2]+(i-x[2])/(x[2]-x[1])*(y[2]-y[1]);
				 if((j-f)*(x[2]-x[1])<0)
				 {   
					 thead=deltriall(thead,t);//ɾȥ���������
				 }
			 }	 
			 
			 else if((i-x[2])*(y[2]-y[1])>0)
			 {
				 thead=deltriall(thead,t);//ɾȥ���������
			 }
			 }
			 
			 else if( (t->n3==m1->num||t->n3==m2->num) &&(t->n2==m1->num||t->n2==m2->num))
			 {   i=node[t->n1].x;
			 j=node[t->n1].y;
			 
			 if (x[1]!=x[2])//�����������
			 {   
				 f=y[2]+(i-x[2])/(x[2]-x[1])*(y[2]-y[1]);
				 if((j-f)*(x[2]-x[1])<0)
				 {
					 thead=deltriall(thead,t);//ɾȥ���������
				 }
			 }	 
			 
			 else if((i-x[2])*(y[2]-y[1])>0)
			 {
				 thead=deltriall(thead,t);//ɾȥ���������
			 }
			 }
		 }
		 
		 
		 
		 
			  }
       

		 printf("shan chu li mian tri\n");
		 //ɾ�������������������
		 p=head->next;//��������Ľڵ�����
		 for(p1=head->next->next;node[p1->num].x!=node[p->num].x || node[p1->num].y!=node[p->num].y;p1=p1->next)
		 {};
		 p1->num=p->num;
          for(m1=p;m1!=p1;m1=m1->next)//�����ڵ�

		 {  
			 m2=m1->next;
			 b1=m1->num;
			 b2=m2->num;

			 x[1]=node[b1].x;
			 y[1]=node[b1].y;
			 float e=0;
			 e=node[b2].x;
			 y[2]=node[b2].y;
			 x[2]=e;


	        	 for (t=thead;t!=NULL;t=t->next)//����������//������һ�У�x[2]��ֵ�仯�ˣ����ˣ�x[2]�ĵ�ַ��t��ͻ��
				 {    							

					 a1=t->n1;
					 a2=t->n2;
					 a3=t->n3;
					 x[1]=node[b1].x;
					 y[1]=node[b1].y;
       				 x[2]=node[b2].x;
					 y[2]=node[b2].y;


					 if( (a1==b1||a1==b2) &&(a2==b1||a2==b2))
						 
					 {   					
						 x[1]=node[b1].x;
						 y[1]=node[b1].y;
						 x[2]=node[b2].x;
						 y[2]=node[b2].y;
						 i=node[a3].x;
						 j=node[a3].y;
						 
						 
						 if (x[1]!=x[2])//�����������
						 {   
							 f=y[2]+(i-x[2])/(x[2]-x[1])*(y[2]-y[1]);
							 if((j-f)*(x[2]-x[1])<0)
							 {
								 thead=deltriall(thead,t);//ɾȥ���������
							 }
						 }	 
						 
						 else if((i-x[2])*(y[2]-y[1])>0)
						 { 
							 thead=deltriall(thead,t);//ɾȥ���������
						 }
					 }
					 else if( (a1==b1||a1==b2) &&(a3==b1||a3==b2))
					 {  
						 i=node[a2].x;
						 j=node[a2].y;
						 
						 
						 if (x[1]!=x[2])//�����������
						 {   
							 f=y[2]+(i-x[2])/(x[2]-x[1])*(y[2]-y[1]);
							 if((j-f)*(x[2]-x[1])<0)
							 {
								 thead=deltriall(thead,t);//ɾȥ���������
							 }
						 }	 
						 
						 else if((i-x[2])*(y[2]-y[1])>0)
						 {
							 thead=deltriall(thead,t);//ɾȥ���������
						 }
					 }
					 else if((a3==b1||a3==b2) &&(a2==b1||a2==b2))
					 {   
					 i=node[t->n1].x;
					 j=node[t->n1].y;
					 
					 
					 if (x[1]!=x[2])//�����������
					 {
						 f=y[2]+(i-x[2])/(x[2]-x[1])*(y[2]-y[1]);
						 
						 if((j-f)*(x[2]-x[1])<0)
						 {
							 thead=deltriall(thead,t);//ɾȥ���������
						 }
					 }	 
					 
					 else if((i-x[2])*(y[2]-y[1])>0)
					 {   
						 thead=deltriall(thead,t);//ɾȥ���������
					 }
					 }

			 
		 }
		 
		 
		 
		 
		 }
return thead;
 }



 void dengfenbianjie(Node(&node)[100000],Node(&mode)[1000][1000],int &m,int &h,float &b,int NN)//����ȱ������α߳�
 {
 int i=0,j=0,n=0; //n���û��涨������
 int m1[1000];
 float m2[1000];//�涨ÿ���ߵȷֵ���Ŀ,��m2ȡ�����õ�m1
 float s;//�����������
 float sb;//�����׼���������
 n=NN;


s=node[m].x*node[1].y-node[1].x*node[m].y;                //���������
for (i=1;i<=m-1;i++)
{s=s+node[i].x*node[i+1].y-node[i+1].x*node[i].y;}
s=s/2;
printf("������������ %f\n",s);

sb=s/n;
printf("����������� %f\n",sb);



b=sqrt((4/1.732)*s/n);
printf("ƽ�������α߳���  %f\n",b);

 for (i=1;i<=m-1;i++)                                        //��m���ߵȷ�
 {m2[i]=d(node[i].x,node[i].y,node[i+1].x,node[i+1].y)/b;
  m1[i]=int(d(node[i].x,node[i].y,node[i+1].x,node[i+1].y)/b);
 }



                       
  for(i=1;i<=m-1;i++)//��m���ߵĵȷֵ�����
  {for(j=1;j<m1[i];j++)
	  {mode[i][j].x=node[i].x+(j)*(node[i+1].x-node[i].x)/m1[i];
	   mode[i][j].y=node[i].y+(j)*(node[i+1].y-node[i].y)/m1[i];
	  }
  }
     for(i=1;i<=m-1;i++)
  {	  for(j=1;j<=m1[i]-1;j++)
      {h++;            
	   node[m+h].x=mode[i][j].x;
	   node[m+h].y=mode[i][j].y;
	  }
  }

 }




void writefile( Node(&node)[100000],int &h,struct tri *thead, FILE *fp2)//д�ļ����ڵ�����ʦ�����
 {    struct tri *t;
      int i=1;
	  fp2=fopen("d:\\tri.txt","w");     //�������������ʼ�������нڵ�.txt
	  fprintf(fp2,"%d\n",h);
	  for(i=1;i<=h;i++)
		  fprintf(fp2,"%d   %f  %f\n",i,node[i].x,node[i].y);
	  fprintf(fp2,"%d\n",count(thead));
	  for(t=thead;t!=NULL;t=t->next)
		  fprintf(fp2,"%d  %d  %d\n",t->n1,t->n2,t->n3);//t->num,

	fclose (fp2);
}


struct tri *bianjietri(struct tri *thead,Node(&node)[100000],int &m,int &h)//���ñ߽�ڵ�����������
{
	struct tri *p,*q,*s,*r;
	int a,b,c,e,f,g,i,k,l;
	int online(Node(&node)[100000],int a,int b,int c);
	int max(struct tri *thead);

		for(p=thead;p!=NULL;p=p->next)
		{ };
		l=max(thead)+1;


	for(p=thead;p->num!=l;p=p->next)//��һ��ͨ���Բ��У�5.9����
	{  
		a=p->n1;b=p->n2;c=p->n3;
	    e=p->t1;f=p->t2;g=p->t3;

		for(i=m+1;i<=m+h;i++)
		{	       
			        if ((online(node,1,2,i)==0)&&(online(node,a,b,i)==1||
						online(node,a,c,i)==1||online (node,b,c,i)==1))//������������α���
					{       
						do                       //��������p��ת������λ��
		                          {k=p->n1;       
		                          p->n1=p->n2;
		                          p->n2=p->n3;
		                          p->n3=k;
	                  	          k=p->t1;       
		                          p->t1=p->t2;
		                          p->t2=p->t3;
		                          p->t3=k;
								  }while(online (node,p->n2,p->n3,i)==0);
								  a=p->n1;b=p->n2;c=p->n3;
					        	if( 
									(online(node,1,2,i+1)==1|| 
									      (online(node,a,b,i+1)==0&&online(node,a,c,i+1)==0&&online (node,b,c,i+1)==0
										  )
								    ) 
						        	||(
									    (online(node,1,2,i+1)==0)&& (online(node,i,b,i+1)==1 )
									  )
								  )	//�����һ���㲻����������α��ϻ����µ�p����i+1
						
							
							
						
								{     
								        e=p->t1;f=p->t2;g=p->t3;

								   for(q=thead;q->next!=NULL;q=q->next)   
								   {};
								  s=(struct tri *)malloc (len);
								  {
									  s->num=max(thead)+1;          //�ڶ�β�����Ҳ��һ��������
									  s->n1=p->n1;
									  s->n2=i;
									  s->n3=p->n3;
									  s->t1=p->num;
									  s->t2=p->t2;
									  s->t3=p->t3;
								  }
								  q->next=s;
								  s->next=NULL;
								  
								  p->n3=i;                        //�����������num ���䣬���¸�ֵ
								  p->t3=s->num;
								  
								
								  
								  //�����ٱ�������
								     if (g!=0)
									 {  
									  for(r=thead;r->num!=g;r=r->next)
									  {};
									  if      (r->t1==p->num) r->t1=s->num;
									  else if (r->t2==p->num) r->t2=s->num;
									  else if (r->t3==p->num) r->t3=s->num;
						
									 }


									 
                        
								}
						else if ((online(node,1,2,i+1)==0)&&(online(node,i,c,i+1)==1))//P���Ҳ�
							{ e=p->t1;f=p->t2;g=p->t3;

								   for(q=thead;q->next!=NULL;q=q->next)   
								   {};
								  s=(struct tri *)malloc (len);
								  {
									  s->num=max(thead)+1;          //�ڶ�β�����Ҳ��һ��������
									  s->n1=p->n1;
									  s->n2=p->n2;
									  s->n3=i;
									  s->t1=p->t1;
									  s->t2=p->t2;
									  s->t3=p->num;
								  }
								  q->next=s;
								  s->next=NULL;
								  
								  p->n2=i;                        //�����������num ���䣬���¸�ֵ
								  p->t1=s->num;
								
								     if (e!=0)
									 {  
									  for(r=thead;r->num!=e;r=r->next)
									  {};
									  if      (r->t1==p->num) r->t1=s->num;
									  else if (r->t2==p->num) r->t2=s->num;
									  else if (r->t3==p->num) r->t3=s->num;

									 }
						}


					}
				
		    
		}
	}
	return thead;
}

int online(Node(&node)[100000],int a,int b,int c)//�ж�c�Ƿ����߶�ab��,����ֵΪ1��ʾ�����߶��ϣ�0��ʾ�㲻������
{
	if(((node[c].x-node[a].x)*(node[b].x-node[c].x))>=0 && ((node[c].y-node[a].y)*(node[b].y-node[c].y))>=0)
	{
		if (node[a].x!=node[b].x)
		{
			if (fabs((node[c].y-node[a].y)*(node[b].x-node[a].x)-(node[c].x-node[a].x)*(node[b].y-node[a].y))<=0.0001)
			      return 1;
		    else  return 0;
		}
    	else if( fabs(node[c].x-node[a].x)<=0.00001)
	     	return 1;
    	else return 0;
	}
	else return 0;

}



void delaunay(Node(&node)[100000],struct tri *thead)
{
   struct tri *d1,*d2,*p;
   int m,n,i,a,b,c,a1,b1,c1,m1,m2,m3;
   float x1,x2,x3,y1,y2,y3,x,y,r;
  
DELAUNAY: for(d1=thead;d1!=NULL;d1=d1->next)
   {   
                x1=node[d1->n1].x; 
                x2=node[d1->n2].x;   
                x3=node[d1->n3].x; 
                y1=node[d1->n1].y; 
                y2=node[d1->n2].y; 
                y3=node[d1->n3].y;   
				x=((y2-y1)*(y3*y3-y1*y1+x3*x3-x1*x1)-(y3-y1)*(y2*y2-y1*y1+x2*x2-x1*x1))/(2*(x3-x1)*(y2-y1)-2*((x2-x1)*(y3-y1)));  
                y=((x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1)-(x3-x1)*(x2*x2-x1*x1+y2*y2-y1*y1))/(2*(y3-y1)*(x2-x1)-2*((y2-y1)*(x3-x1)));  
                r=sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y));


	   	 m1=d1->t1;m2=d1->t2;m3=d1->t3;	
		 
	   if (m1!=0)
	   { 
		   for(d2=thead;d2->num!=m1;d2=d2->next)
		   {};

		  do                       //��d1��������ת���������λ��
		  {n=d1->n1;       
		  d1->n1=d1->n2;
		  d1->n2=d1->n3;
		  d1->n3=n;
		  m=d1->t1;       
		  d1->t1=d1->t2;
		  d1->t2=d1->t3;
		  d1->t3=m;
		  }while(d1->n2==d2->n1 || d1->n2==d2->n2  ||  d1->n2==d2->n3);

		  do                       //��d2��������ת���������λ��
		  {n=d2->n1;       
		  d2->n1=d2->n2;
		  d2->n2=d2->n3;
		  d2->n3=n;
		  m=d2->t1;       
		  d2->t1=d2->t2;
		  d2->t2=d2->t3;
		  d2->t3=m;
		  }while(d2->n2==d1->n1 || d2->n2==d1->n2  ||  d2->n2==d1->n3);

		  a=d1->n1; b=d1->n2; c=d1->n3;
		  a1=d1->t1;b1=d1->t2;c1=d1->t3;
              
		  i=d2->n2;
		   if ((((node[i].x-x)*(node[i].x-x)+(node[i].y-y)*(node[i].y-y)-r*r))<0 && 
			   fabs((node[i].x-x)*(node[i].x-x)+(node[i].y-y)*(node[i].y-y)-r*r)>0.0002)//���d2ָ�����������d1���Բ��
		   {  
			   d1->n1=d2->n2;
			   d1->n2=a;
			   d1->n3=b;

			   d2->n1=b;
			   d2->n3=d2->n2;
			   d2->n2=c;

			   d1->t1=d2->t2;
			   d1->t2=a1;
			   
			   d2->t2=d2->t1;
			   d2->t1=b1;    //������������n1,n2,n3,t1,t2,t3���跭��֮�����ֵ

               m=d1->t1;
			   if(m!=0)
			   {
				   for (p=thead;p->num!=m;p=p->next)
				   {};
				   if (p->t1==d2->num)      p->t1=d1->num;
				   else if (p->t2==d2->num) p->t2=d1->num;
				   else if (p->t3==d2->num) p->t3=d1->num;

			   }

			   m=d2->t1;
			   if(m!=0)
			   {
				   for (p=thead;p->num!=m;p=p->next)
				   {};
				   if (p->t1==d1->num)      p->t1=d2->num;
				   else if (p->t2==d1->num) p->t2=d2->num;
				   else if (p->t3==d1->num) p->t3=d2->num;

			   }


           goto DELAUNAY;
		   }


	   }

	    if (m2!=0)//��һ�����ˣ���ת֮��ԭ���������Ѿ��仯��
	   {
		   for(d2=thead;d2->num!=m2;d2=d2->next)
		   {};

		  do                       //��d1��������ת���������λ��
		  {n=d1->n1;       
		  d1->n1=d1->n2;
		  d1->n2=d1->n3;
		  d1->n3=n;
		  m=d1->t1;       
		  d1->t1=d1->t2;
		  d1->t2=d1->t3;
		  d1->t3=m;
		  }while(d1->n2==d2->n1 || d1->n2==d2->n2  ||  d1->n2==d2->n3);

		  do                       //��d2��������ת���������λ��
		  {n=d2->n1;       
		  d2->n1=d2->n2;
		  d2->n2=d2->n3;
		  d2->n3=n;
		  m=d2->t1;       
		  d2->t1=d2->t2;
		  d2->t2=d2->t3;
		  d2->t3=m;
		  }while(d2->n2==d1->n1 || d2->n2==d1->n2  ||  d2->n2==d1->n3);
           
		  a=d1->n1; b=d1->n2; c=d1->n3;
		  a1=d1->t1;b1=d1->t2;c1=d1->t3;
		  i=d2->n2;
		   if ((((node[i].x-x)*(node[i].x-x)+(node[i].y-y)*(node[i].y-y)-r*r))<0 && 
			   fabs((node[i].x-x)*(node[i].x-x)+(node[i].y-y)*(node[i].y-y)-r*r)>0.0002)
		   {    
			   d1->n1=a; d1->n2=b; d1->n3=c;
		       d1->t1=a1;d1->t2=b1;d1->t3=c1;

			   d1->n1=d2->n2;
			   d1->n2=a;
			   d1->n3=b;

			   d2->n1=b;
			   d2->n3=d2->n2;
			   d2->n2=c;

			   d1->t1=d2->t2;
			   d1->t2=a1;
			   
			   d2->t2=d2->t1;
			   d2->t1=b1;    //������������n1,n2,n3,t1,t2,t3���跭��֮�����ֵ

               m=d1->t1;
			   if(m!=0)
			   {
				   for (p=thead;p->num!=m;p=p->next)
				   {};
				   if (p->t1==d2->num)      p->t1=d1->num;
				   else if (p->t2==d2->num) p->t2=d1->num;
				   else if (p->t3==d2->num) p->t3=d1->num;

			   }

			   m=d2->t1;
			   if(m!=0)
			   {
				   for (p=thead;p->num!=m;p=p->next)
				   {};
				   if (p->t1==d1->num)      p->t1=d2->num;
				   else if (p->t2==d1->num) p->t2=d2->num;
				   else if (p->t3==d1->num) p->t3=d2->num;

			   }

             goto DELAUNAY;

		   }


	   }

	   if (m3!=0)//���Ӧ���Ƕ���
	   {
		   for(d2=thead;d2->num!=m3;d2=d2->next)
		   {};

		  do                       //��d1��������ת���������λ��
		  {n=d1->n1;       
		  d1->n1=d1->n2;
		  d1->n2=d1->n3;
		  d1->n3=n;
		  m=d1->t1;       
		  d1->t1=d1->t2;
		  d1->t2=d1->t3;
		  d1->t3=m;
		  }while(d1->n2==d2->n1 || d1->n2==d2->n2  ||  d1->n2==d2->n3);

		  do                       //��d2��������ת���������λ��
		  {n=d2->n1;       
		  d2->n1=d2->n2;
		  d2->n2=d2->n3;
		  d2->n3=n;
		  m=d2->t1;       
		  d2->t1=d2->t2;
		  d2->t2=d2->t3;
		  d2->t3=m;
		  }while(d2->n2==d1->n1 || d2->n2==d1->n2  ||  d2->n2==d1->n3);
           
		  a=d1->n1; b=d1->n2; c=d1->n3;
		  a1=d1->t1;b1=d1->t2;c1=d1->t3;
		  i=d2->n2;
		   if ((((node[i].x-x)*(node[i].x-x)+(node[i].y-y)*(node[i].y-y)-r*r))<0 && 
			   fabs((node[i].x-x)*(node[i].x-x)+(node[i].y-y)*(node[i].y-y)-r*r)>0.0002)
		   {   
			   d1->n1=a; d1->n2=b; d1->n3=c;
		       d1->t1=a1;d1->t2=b1;d1->t3=c1;

			   d1->n1=d2->n2;
			   d1->n2=a;
			   d1->n3=b;

			   d2->n1=b;
			   d2->n3=d2->n2;
			   d2->n2=c;

			   d1->t1=d2->t2;
			   d1->t2=a1;
			   
			   d2->t2=d2->t1;
			   d2->t1=b1;    //������������n1,n2,n3,t1,t2,t3���跭��֮�����ֵ

               m=d1->t1;
			   if(m!=0)
			   {
				   for (p=thead;p->num!=m;p=p->next)
				   {};
				   if (p->t1==d2->num)      p->t1=d1->num;
				   else if (p->t2==d2->num) p->t2=d1->num;
				   else if (p->t3==d2->num) p->t3=d1->num;

			   }

			   m=d2->t1;
			   if(m!=0)
			   {
				   for (p=thead;p->num!=m;p=p->next)
				   {};
				   if (p->t1==d1->num)      p->t1=d2->num;
				   else if (p->t2==d1->num) p->t2=d2->num;
				   else if (p->t3==d1->num) p->t3=d2->num;

			   }

            goto DELAUNAY;

		   }


	   }



   }
}


 






void insert(struct tri *thead,Node(&node)[100000],float b,int &h)//��ʽ���ʼ�����ڲ�����㣬���ɵȱ�������ֱ�����������ζ��㹻С.b�������α�׼�߳�
{   struct tri *p,*p0,*p1,*o;
    float x[4],y[4],f[4];
	float i3,i4,j3,j4,mm;//mm��ʾ�д���б��
	int i=0,j=0;
	int k;
	int satisfy(Node(&node)[100000],struct tri *p,float &b);
	int l1,l2;//�洢intriornot�ĺ���ֵ
	int n;//��¼�����ɵ����������ε��ٱ�
    int ntri; //��¼����������������
	int liu=0;
for(liu=0;liu<=0;liu++)
  {	for(p=thead;p!=NULL;p=p->next)

	{    ntri=max(thead);
		if(satisfy(node,p,b)==0)
		{   	
			    x[1]=node[p->n1].x;
		        x[2]=node[p->n2].x;
				x[3]=node[p->n3].x;
				y[1]=node[p->n1].y;
				y[2]=node[p->n2].y;
				y[3]=node[p->n3].y;
				x[4]=node[p->n1].x;
				y[4]=node[p->n1].y;
           
			if( fabs(d(x[1],y[1],x[2],y[2])-b)<=0.1||
				fabs(d(x[1],y[1],x[3],y[3])-b)<=0.1||
				fabs(d(x[3],y[3],x[2],y[2])-b)<=0.1 )
			{  
				 do                       //��������p��ת������λ��,pͬn2,n3�γɵȱ�������
				 {k=p->n1;       
				 p->n1=p->n2;
				 p->n2=p->n3;
				 p->n3=k;
				 k=p->t1;       
				 p->t1=p->t2;
				 p->t2=p->t3;
				 p->t3=k;
				 x[1]=node[p->n1].x;
				 x[2]=node[p->n2].x;
				 x[3]=node[p->n3].x;
				 y[1]=node[p->n1].y;
				 y[2]=node[p->n2].y;
				 y[3]=node[p->n3].y;
				 x[4]=x[1];
				 y[4]=y[1];

				 }while(fabs(d(x[1],y[1],x[2],y[2])-b)>=0.1);//�˵ز���
				 
				 if (fabs(y[2]-y[1])>0.0001)
				 {
					 mm=-(x[2]-x[1])/(y[2]-y[1]);
					 i3=0.5*(x[1]+x[2])+sqrt(0.75*b*b/(1+mm*mm));
					 i4=0.5*(x[1]+x[2])-sqrt(0.75*b*b/(1+mm*mm));
					 j3=mm*(i3-0.5*x[1]-0.5*x[2])+0.5*(y[1]+y[2]);
					 j4=mm*(i4-0.5*x[1]-0.5*x[2])+0.5*(y[1]+y[2]);
				 }
				 else if(fabs(y[2]-y[1])<=0.0001)
				 {
					 i3=0.5*(x[1]+x[2]);
					 i4=0.5*(x[1]+x[2]);
					 j3=y[1]+sqrt(3)*0.5*b;
					 j4=y[1]-sqrt(3)*0.5*b;
				 } 

				 for(j=0,i=1;i<=3;i++)                           //�ҵ�������a������nodenum
				 {   
					 if (x[i]!=x[i+1])
						 
					 {   
						 f[i]=y[i]+(i3-x[i])/(x[i+1]-x[i])*(y[i+1]-y[i]);
						 if((j3-f[i])*(x[i+1]-x[i])>0)
						 {
							 j++;
						 }
					 }
					 else if((i3-x[i+1])*(y[i+1]-y[i])<0)
					 {j++;}
				 }
					 if (j==3)
					 { 
						l1=1;;
					 }
					 
					 else
					 {   
						 l1=0;
						 
					 }
				 
					 for(j=0,i=1;i<=3;i++)                           //�ҵ�������a������nodenum
					 {   
						 if (x[i]!=x[i+1])
							 
						 {   
							 f[i]=y[i]+(i4-x[i])/(x[i+1]-x[i])*(y[i+1]-y[i]);
							 if((j4-f[i])*(x[i+1]-x[i])>0)
							 {
								 j++;
							 }
						 }
					 
						 else if((i4-x[i+1])*(y[i+1]-y[i])<0)
						 {j++;}
					 }	 
						 if (j==3)
						 { 
							 l2=1;
						 }
						 
						 else
						 {   
							 l2=0;
							 
						 }
					 
  
	
				   	

					 if(l1==1 || l2==1)
					 {  
						 if(l1==1)
						 {   
							 h=h+1;
							 node[h].x=i3;
							 node[h].y=j3;

						 }
						 else if(l2==1)
						 {
							 h=h+1;
							 node[h].x=i4;
							 node[h].y=j4;
							 
						 }
						 
						 for(o=thead ;o->next!=NULL;o=o->next)//��ָ��dָ���β������
						 {};   
						 p0=(struct tri *)malloc (len);
						 {      
							 p0->num=ntri+1;          //�ڶ�β�����һ��������
							 p0->n1=h;
							 p0->n2=p->n2;
							 p0->n3=p->n3;
							 p0->t1=p->num;
							 p0->t2=p->t2;
							 p0->t3=ntri+2;
							 o->next=p0;
							 p0->next=NULL;
						 
						 }     
						 p1=(struct tri *)malloc (len);
						 
						 {
							 p1->num=ntri+2;          //�ڶ�β����ڶ���������
							 p1->n1=p->n3;
							 p1->n2=p->n1;
							 p1->n3=h;
							 p1->t1=p->t3;
							 p1->t2=p->num;
							 p1->t3=ntri+1;
							 p0->next=p1;
							 p1->next=NULL;
						 }
						 
						 {
							 p->n3=h;         //�滻ԭ��������
							 p->t2=ntri+1;
							 p->t3=ntri+2;
						 }
						 
						 n=p0->t2;//�ҵ���һ�������εڶ������ڱ�t2����Ӧ����������ţ����ı���t2��ֵ������Ϊ��ʹԭ��������������Χ���ĸ�������t1,t2,t3��ȷ��
						 if (n!=0)
						 {
							 for(o=thead;o->num!=n;o=o->next)
							 {};
							 if (o->t1==p->num)          o->t1=(ntri+1);
							 else if (o->t2==p->num)     o->t2=(ntri+1);
							 else if (o->t3==p->num)     o->t3=(ntri+1);
						 }
						 
						 n=p1->t1;//�ҵ���һ�������εڶ������ڱ�t2����Ӧ����������ţ����ı���t2��ֵ������Ϊ��ʹԭ��������������Χ���ĸ�������t1,t2,t3��ȷ��
						 if (n!=0)
						 {
							 for(o=thead;o->num!=n;o=o->next)
							 {};
							 if (o->t1==p->num)          o->t1=(ntri+2);
							 else if (o->t2==p->num)     o->t2=(ntri+2);
							 else if (o->t3==p->num)     o->t3=(ntri+2);

						 }
						 
					 }




					 
					 
					 
			}
			
			
			
		}
		
	}
delaunay(node,thead);
}

}

 int satisfy(Node(&node)[100000],struct tri *p,float &b)//�ж�������p�Ƿ����������߶��ӽ��ڱ߳�b��Ҫ�����㷵��1�������㷵��0
  {
	  float f=0;
	  float x[3],y[3],r[3];
	        x[1]=node[p->n1].x;
			x[2]=node[p->n2].x;
			x[3]=node[p->n3].x;
			y[1]=node[p->n1].y;
			y[2]=node[p->n2].y;
			y[3]=node[p->n3].y;
			r[1]=sqrt((x[1]-x[2])*(x[1]-x[2])+(y[1]-y[2])*(y[1]-y[2]));
			r[2]=sqrt((x[3]-x[2])*(x[3]-x[2])+(y[3]-y[2])*(y[3]-y[2]));
			r[3]=sqrt((x[1]-x[3])*(x[1]-x[3])+(y[1]-y[3])*(y[1]-y[3]));
			f=sqrt((r[1]-b)*(r[1]-b)+(r[2]-b)*(r[2]-b)+(r[3]-b)*(r[3]-b));
			if(f<0.5)//����ǳ��ؼ����պ�����ܷ�����ȫ������
				return 1;
			else return 0;
  }

  void sinsert (struct tri *thead,Node(&node)[100000],float &b,int &h)//2.������ܴ���������ڲ���㣬����ʣ��������
 {   
	 struct tri *p;
	 p=thead;
	 float S=0,	SS;
	 float x[4];
	 float y[4];
	 int k;//i���ڼ�����j���ڱ����ж��ٸ������β��ϸ�,k������ת������������λ��
	 float l[4];//�����������ߵĳ���
	 float r[3];//�����������ǵĶ���
	 float i,j;//��ʾ���������ԲԲ�ĺ�������
	 int trinum[1000],ntri,nodenum;
	 int cc;
	 for(cc=0;cc<10;cc++)
	 

     {   
    	 for (p=thead;p!=NULL;p=p->next)
		 {
		     x[1]=node[p->n1].x;//x[1],x[2],x[3]��ʾ�����ζ������꣬x[4]�����ԲԲ������
			 x[2]=node[p->n2].x;
			 x[3]=node[p->n3].x;
			 y[1]=node[p->n1].y;
			 y[2]=node[p->n2].y;
			 y[3]=node[p->n3].y;
			 S = s(x,y);
			 SS=S/(sqrt(3)/4*b*b);


			   if (SS>1.5 && SS<3)
			   {   

				   
				 
				   if( fabs(d(x[1],y[1],x[2],y[2])/b-1)<=0.2||
					   (fabs(d(x[1],y[1],x[3],y[3])/b-1))<=0.2||
					   fabs(d(x[3],y[3],x[2],y[2])/b-1)<=0.2 )
				   {  
					 do                       //��������p��ת������λ��,pͬn2,n3�γɵȱ�������
					
					 { thead=xuanzhuan(thead, p->num) ;
					 x[1]=node[p->n1].x;
					 x[2]=node[p->n2].x;
					 x[3]=node[p->n3].x;
					 y[1]=node[p->n1].y;
					 y[2]=node[p->n2].y;
					 y[3]=node[p->n3].y;

					 
					 }while(fabs(d(x[1],y[1],x[2],y[2])/b-1)>=0.2);
					 l[1]=d(x[1],y[1],x[2],y[2]);
					 l[2]=d(x[3],y[3],x[2],y[2]);
					 l[3]=d(x[1],y[1],x[3],y[3]);
					 r[1]=acos((l[1]*l[1]+l[3]*l[3]-l[2]*l[2])/(2*l[1]*l[3]));
					 r[2]=acos((l[1]*l[1]+l[2]*l[2]-l[3]*l[3])/(2*l[1]*l[2]));
					 
                     trinum[1]=p->num;
					   if (fabs(3.1415926/3-r[1])>fabs(3.1415926/3-r[2]))
					   { h=h+1;
						 trinum[2]=p->t2;
						 node[h].x=x[2]-b/l[3]*(x[2]-x[3]);
						 node[h].y=y[2]-b/l[3]*(y[2]-y[3]);
					   }
					  else 
					  {   h=h+1;
						  trinum[2]=p->t3;
						  node[h].x=x[1]-b/l[3]*(x[1]-x[3]);
					      node[h].y=y[1]-b/l[3]*(y[1]-y[3]);
					  }
					     nodenum=h;
					 
					      if(trinum[2]!=0)
						  {
						      
							  newd(thead,trinum,ntri,nodenum);

							  
						      delaunay(node,thead);
						 
						  }	 
				 	     
						 
						 }
			  
			   }
		 
		 
		 }
		  
		     for (p=thead;p!=NULL;p=p->next)
			 {
			 
				 x[1]=node[p->n1].x;//x[1],x[2],x[3]��ʾ�����ζ������꣬x[4]�����ԲԲ������
				 x[2]=node[p->n2].x;
				 x[3]=node[p->n3].x;
				 y[1]=node[p->n1].y;
				 y[2]=node[p->n2].y;
				 y[3]=node[p->n3].y;
                 i=(x[1]+x[2]+x[3])/3;
				 j=(y[1]+y[2]+y[3])/3;
				 S = s(x,y);
				 SS=S/(sqrt(3)/4*b*b);
			      if (SS>=3)
				  {   
					  
					  
				      trinum[1]=p->num;
					  h=h+1;
					  node[h].x=i;
					  node[h].y=j;
					  nodenum=h;
					  ntri=max(thead);
					  newi(thead,trinum,ntri,nodenum);
					  
					  delaunay(node,thead);
					  
				  }
			 }
	 }
 }
 float  mintri(struct tri *thead,Node(&node)[100000],int i)//����������̱߳�
 {
	 struct tri *p1;
	 float a1,a2,a3;
	 for(p1=thead;p1->num!=i;p1=p1->next)
	 {};
     a1=d(node[p1->n1].x,node[p1->n1].y,node[p1->n2].x,node[p1->n2].y);
	 a2=d(node[p1->n1].x,node[p1->n1].y,node[p1->n3].x,node[p1->n3].y);
	 a3=d(node[p1->n3].x,node[p1->n3].y,node[p1->n2].x,node[p1->n2].y);
	 if (a1>=a2)
	 {
		 a1=a2;
	 }
	 if (a1>=a3)
	 {
		 a1=a3;
	 }
	 return a1;


 }
 float  maxtri(struct tri *thead,Node(&node)[100000],int i)//����������߳�
 {
	 struct tri *p1;
	 float a1,a2,a3;
	 for(p1=thead;p1->num!=i;p1=p1->next)
	 {};
     a1=d(node[p1->n1].x,node[p1->n1].y,node[p1->n2].x,node[p1->n2].y);
	 a2=d(node[p1->n1].x,node[p1->n1].y,node[p1->n3].x,node[p1->n3].y);
	 a3=d(node[p1->n3].x,node[p1->n3].y,node[p1->n2].x,node[p1->n2].y);
	 if (a1<=a2)
	 {
		 a1=a2;
	 }
	 if (a1<=a3)
	 {
		 a1=a3;
	 }
	 return a1;
	 
	 
 }
 
void youhua(struct tri *thead,Node(&node)[100000],float &b)//ɾ��̫С�������Σ��������ڵ�ϲ���һ������������������
{    
	 struct tri *p,*p1,*p2,*q,*f;
	 int i,j;
	 
	 for(p1=thead;p1!=NULL;p1=p1->next)
	 {
		 
		 if(mintri(thead,node,p1->num)<=(0.3*b))
		 {   printf("p1->num=%d\n",p1->num);
			 do
			 {
				 thead=xuanzhuan(thead,p1->num);
				 
			 }while (d(node[p1->n3].x,node[p1->n3].y,node[p1->n2].x,node[p1->n2].y)>=(0.3*b));//��������p1��ת������λ��

			 if(p1->t2!=0)
			 {
				 for(p2=thead;p2->num!=p1->t2;p2=p2->next)
				 {};
			 printf("p2->num=%d\n",p2->num);
				 do
				 {
					 thead=xuanzhuan(thead, p2->num);
				 }
				 while (d(node[p2->n3].x,node[p2->n3].y,node[p2->n2].x,node[p2->n2].y)>=(0.3*b));

				 node[p1->n2].x=0.5*(node[p1->n2].x+node[p1->n3].x);
				 node[p1->n2].y=0.5*(node[p1->n2].y+node[p1->n3].y);
				 node[p1->n3].x=0;
			     node[p1->n3].y=0;
				
				 //����������6
				 if(p1->t1!=0)
				 {
					 for(q=thead;q->num!=p1->t1;q=q->next)
					 {};
					 do
					 {
						 thead=xuanzhuan(thead,q->num);

					 }while (q->t1!=p1->num);
					 q->t1=p1->t3;

				 }
				 //����������5
				 if(p1->t3!=0)
				 {
				 
					 for(q=thead;q->num!=p1->t3;q=q->next)
					 {};
					 do
					 {
						 thead=xuanzhuan(thead,q->num);
						 
					 }while (q->t1!=p1->num);
					 q->t1=p1->t1;
					 q->n2=p1->n2;

				 }
				 //����������2,ͬ5��·��ͨ������t2
				 if(p2->t3!=0)
				 {
					 
					 for(q=thead;q->num!=p2->t3;q=q->next)
					 {};
					 do
					 {
						 thead=xuanzhuan(thead,q->num);
						 
					 }while (q->t1!=p2->num);
					 q->t1=p2->t1;
	
					 
				 }
				 //����������3
				 if(p2->t1!=0)
				 {
					 for(q=thead;q->num!=p2->t1;q=q->next)
					 {};
					 do
					 {
						 thead=xuanzhuan(thead,q->num);
						 
					 }while (q->t1!=p2->num);
					 q->t1=p2->t3;
					 q->n1=p1->n2;
					 
				 }
				 //����������4�Լ���һ�������λ�õ�������
				 i=p1->n3;
				 j=p1->n2;
				 thead=deltri (thead,p1->num);
				 thead=deltri (thead,p2->num);
				 for(p=thead;p!=NULL;p=p->next)
				 {
					 if(p->n1==i||p->n2==i||p->n3==i)
					 {  
						 do
						 {
							 thead=xuanzhuan(thead,p->num);
							 
					  }while(p->n1!=i);
						 p->n1=j;

					 }
				 }




			 
			 }
			 
		 }

	 }

 delaunay(node,thead);
 }

  void lapace(struct tri *thead,Node(&node)[100000],int &hh,int &h)//lapace �Ż�����
  {
  
	  int i=1,j=1,k=1,n=1;
	  float x[1000000],y[1000000],XX=0,YY=0;//x�洢����������ڵ������֮�ͣ�y�洢������֮�ͣ�XX,YY���ۼ�֮��
	  float w=1;//w���ɳ����ӣ��������ѡ�����
	  struct tri*p;
	  i=(hh+1);
	  int t=1;
    //  for(t=1;t<=100;t++)	  
	  {for(i=hh+1;i<=h;i++)
	  {
		  if( (node[i].x*node[i].x+node[i].y*node[i].y)!=0)
		  {    j=1;
		       XX=0;
			   YY=0;
		      for(p=thead;p!=NULL;p=p->next)
			  {
				  if(p->n1==i||p->n2==i||p->n3==i)
				  {  
					  do
					  {
						  thead=xuanzhuan(thead,p->num);
						  
					  }while(p->n1!=i);
					  x[j]=node[p->n2].x+node[p->n3].x;
					  y[j]=node[p->n2].y+node[p->n3].y;
					  j=j+1;
				  }  
				        
			  }
			  n=j-1;
			  for(k=1;k<=n;k++)
			  {
				  XX=XX+x[k];
				  YY=YY+y[k];
			  }
			  node[i].x=node[i].x+w*(XX/(2*n)-node[i].x);
			  node[i].y=node[i].y+w*(YY/(2*n)-node[i].y);

		  }
		  
	  

	  }
	  }

  }
 



  void WRitefile( Node(&node)[100000],int &h,struct tri *thead, FILE *fp3)//д�ļ����ڵ�����ʦ�����
  {    struct tri *t;
  int i=1;
  fp3=fopen("d:\\numtri.txt","w");     //�������������ʼ�������нڵ�.txt
  fprintf(fp3,"%d\n",h);
  for(i=1;i<=h;i++)
	  fprintf(fp3,"%d   %f  %f\n",i,node[i].x,node[i].y);
  fprintf(fp3,"%d\n",count(thead));
  for(t=thead;t!=NULL;t=t->next)
	  fprintf(fp3,"%d %d  %d  %d\n",t->num,t->n1,t->n2,t->n3);//t->num,
  
  fclose (fp3);
}

  void deinsert (struct tri *thead,Node(&node)[100000],float &b,int &h)//�����������ģ��
  {   
	  struct tri *p;
	  p=thead;
	  float S=0,	SS;
	  float x[4];
	  float y[4];
	  int k;//i���ڼ�����j���ڱ����ж��ٸ������β��ϸ�,k������ת������������λ��
	  float l[3];//�����������ߵĳ���
	  float r[3];//�����������ǵĶ���
	  float i,j;//��ʾ���������ԲԲ�ĺ�������
	  int trinum[1000],ntri,nodenum;
	  int cc;
//	  for(cc=0;cc<10;cc++)
		  
		  
	  {   
		  for (p=thead;p!=NULL;p=p->next)
		  {
			  if(maxtri(thead,node,p->num)/mintri(thead,node,p->num)>=3.5||mintri(thead,node,p->num)>=(1.5*b)||maxtri(thead,node,p->num)>=2*b)
			  {		 
				  x[1]=node[p->n1].x;//x[1],x[2],x[3]��ʾ�����ζ������꣬x[4]�����ԲԲ������
				  x[2]=node[p->n2].x;
				  x[3]=node[p->n3].x;
				  y[1]=node[p->n1].y;
				  y[2]=node[p->n2].y;
				  y[3]=node[p->n3].y;
				  i=(x[1]+x[2]+x[3])/3;
				  j=(y[1]+y[2]+y[3])/3;
				  trinum[1]=p->num;
				  h=h+1;
				  node[h].x=i;
				  node[h].y=j;
				  nodenum=h;
				  ntri=max(thead);
				  newi(thead,trinum,ntri,nodenum);
				  
				  delaunay(node,thead);
			  }
			  
			  
		  }
		  
		  
	  }
	}


  void qualitytri(Node(&node)[100000],struct tri *thead)//ͳ��������������ָ��
  {   
	  FILE *fp;
	  float x[4],y[4],l[4],r[4];
	  fp=fopen("d:\\quality.txt","w");
	  struct tri *p;
	  float a;
	  for(p=thead;p!=NULL;p=p->next)
	  {
		  
		  x[1]=node[p->n1].x;
		  x[2]=node[p->n2].x;
		  x[3]=node[p->n3].x;
		  y[1]=node[p->n1].y;
		  y[2]=node[p->n2].y;
		  y[3]=node[p->n3].y;
	      l[1]=d(x[1],y[1],x[2],y[2]);
		  l[2]=d(x[3],y[3],x[2],y[2]);
		  l[3]=d(x[1],y[1],x[3],y[3]);
		  r[1]=acos((l[1]*l[1]+l[3]*l[3]-l[2]*l[2])/(2*l[1]*l[3]));
		  r[2]=acos((l[1]*l[1]+l[2]*l[2]-l[3]*l[3])/(2*l[1]*l[2]));
		  r[3]=3.1415926-r[1]-r[2];
		  if(r[1]>=r[2])
                r[1]=r[2];
		  if(r[1]>=r[3])
			    r[1]=r[3];
		  r[1]=r[1]*3/3.1415926;
		  a= maxtri(thead,node,p->num)/mintri(thead,node,p->num);
		  fprintf(fp,"%f  %f\n",a,r[1]);
		  
	  }
	  
	  
  }
