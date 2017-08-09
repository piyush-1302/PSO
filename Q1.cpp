#include<bits/stdc++.h>
using namespace std;
#define ll long long int
#define er 1e-200
#define Rand(low,upper) (((double)rand()*(upper-low)/(double)RAND_MAX)+low)
#define F(x) (x)*(x)
double particles[1007][35],velocity[1007][35],pbest[1007][35],vel[1007][35],gbestar[35],x[1007][35],pbtupd[1007];
double c1,c2,vmax;
ll no_particles,dimensions,lower_limit,upper_limit,gbest_index,max_iter;


void initialise()
{
    for(ll i=0;i<no_particles;i++)
    {
        for(ll j=0;j<dimensions;j++)
        {
            particles[i][j]=Rand(lower_limit,upper_limit);
            velocity[i][j]=0;
        }
    }
}


void firstfitness()
{
    double sum;
    for(ll i=0;i<no_particles;i++)
    {
        sum=0;
        for(ll j=0;j<dimensions;j++)
        {
            sum+=F(particles[i][j]);
        }
        pbtupd[i]=sum;
    }
}

void pbestfinding()
{
    //This function is used for finding pbest/lbest in which RING TOPOLOGY has been used.
   for(ll i=0;i<no_particles;i++)
   {
        if(i==0)
        {
            if(pbtupd[i]<=pbtupd[i+1] && pbtupd[i]<=pbtupd[no_particles-1])
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[i][j];
            }
            else if(pbtupd[i+1]<pbtupd[i] && pbtupd[i+1]<pbtupd[no_particles-1])
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[i+1][j];
            }
            else
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[no_particles-1][j];
            }
        }
        else if(i==no_particles-1)
        {
            if(pbtupd[i]<=pbtupd[0] && pbtupd[i]<=pbtupd[no_particles-2])
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[i][j];
            }
            else if(pbtupd[0]<pbtupd[i] && pbtupd[0]<pbtupd[no_particles-2])
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[0][j];
            }
            else
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[no_particles-2][j];
            }
        }
        else
        {
            if(pbtupd[i]<=pbtupd[i+1] && pbtupd[i]<=pbtupd[i-1])
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[i][j];
            }
            else if(pbtupd[i+1]<pbtupd[i] && pbtupd[i+1]<pbtupd[i-2])
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[i+1][j];
            }
            else
            {
                for(ll j=0;j<dimensions;j++)
                    pbest[i][j]=particles[i-1][j];
            }
        }
   }
}
double firstglobal(double gbest)
{
    //This function is used for finding gbest in which STAR TOPOLOGY has been used.
    ll index;
    for(ll i=0;i<no_particles;i++)
    {
        if(gbest>pbtupd[i])
        {
            gbest=pbtupd[i];
            index=i;
        }
    }
    for(ll j=0;j<dimensions;j++)
        gbestar[j]=particles[index][j];
    return gbest;
}


void update(double w)
{
    for(ll i=0;i<no_particles;i++)
    {
        for(ll j=0;j<dimensions;j++)
        {
            velocity[i][j]=w*velocity[i][j]+c1*Rand(0,1)*(pbest[i][j]-particles[i][j])+c2*Rand(0,1)*(gbestar[j]-particles[i][j]);
            if(velocity[i][j]<-vmax)
                velocity[i][j]=-vmax;
            else if(velocity[i][j]>vmax)
                velocity[i][j]=vmax;

            particles[i][j]=particles[i][j]+velocity[i][j];
            if(particles[i][j]<lower_limit)
                particles[i][j]=lower_limit;
            else if(particles[i][j]>upper_limit)
                particles[i][j]=upper_limit;
        }
    }
}


void fitness()
{
    double sum;
    for(ll i=0;i<no_particles;i++)
    {
        sum=0;
        for(ll j=0;j<dimensions;j++)
        {
            sum+=F(particles[i][j]);
        }
        if(sum<pbtupd[i])
        {
            pbtupd[i]=sum;
        }
    }
}


double better(double gbest)
{
    //This function is used for finding gbest in which STAR TOPOLOGY has been used.
    for(ll i=0;i<no_particles;i++)
    {
        if(gbest>pbtupd[i])
        {
            gbest=pbtupd[i];
            for(ll j=0;j<dimensions;j++)
                gbestar[j]=pbest[i][j];
        }
    }
    return gbest;
}


int main()
{
    ll iter;
    double gbest=100000.0,w,z,zk;
    cin>>no_particles>>dimensions>>lower_limit>>upper_limit>>max_iter>>c1>>c2;//>>vmax;
    initialise();
    firstfitness();
    pbestfinding();
    gbest=firstglobal(gbest);
    iter=1;
    vmax=0.1*(upper_limit-lower_limit);
    zk=0.6;
    while(iter<=max_iter && gbest>=er)
    {
        z=zk*4.0*(1-zk);
        zk=z;
        w=0.5*Rand(0,1)/2+0.5*z;
        update(w);
        fitness();
        pbestfinding();
        gbest=better(gbest);
        iter++;
    }
    cout<<"MINIMUM VALUE : "<<gbest<<"\n";
    return 0;
}
