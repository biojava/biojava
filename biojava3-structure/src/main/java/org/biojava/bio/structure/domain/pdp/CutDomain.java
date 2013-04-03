package org.biojava.bio.structure.domain.pdp;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;


public class CutDomain {

	int ndom;

	List<Domain> domains;

	public static boolean verbose = true;
	
	int[][] dist;
	Atom[] ca;
	
	public CutDomain(Atom[]ca, PDPDistanceMatrix pdpMatrix){
		dist = pdpMatrix.getDist();
		this.ca = ca;

		ndom = 0;
		
		domains = new ArrayList<Domain>();
	
	}


	public  void cutDomain(Domain dom, CutSites cut_sites, PDPDistanceMatrix pdpMatrix){

		if ( verbose )
		System.out.println("  B ... beginning of cutDomain " +dom + " cutsites: " + cut_sites );
		
		/* recursive function to cut input domain into two domains */
		
		
		int i,site;
		
		Domain dom1 = new Domain();
		Domain dom2 = new Domain();
		
		CutValues val = new CutValues();
		val.s_min = 100;
		val.site2 = 0;
		val.first_cut = true;
		
			
		Cut cut = new Cut();

		site = cut.cut(ca,dom,val, dist, pdpMatrix);
		if ( verbose )
		System.out.println("  S ... site " + dom + " : site: " + site + " val : " + val);
		
		if(site<0) {
			
			/* function cut makes a decision where to cut , returns -1 if no cut */
			//memcpy(&domains[ndom],&dom,sizeof(struct Domain));
			domains.add(dom);
			dom.score = val.s_min;
			//dom = domains[ndom];
			ndom++;
			return;
		}
		
		if(verbose) 
			System.out.println(String.format("   C ... Cutting at position(s): %d %d %f\n",site,val.site2,dom.score));
		
		cut_sites.cut_sites[cut_sites.ncuts++] = site;
		
		/* create new domains: dom1 and dom2*/
		dom1.size = 0;
		dom1.nseg = 0;
		dom2.size = 0;
		dom2.nseg = 0;
		if(val.site2==0) { /* single cut*/
			for(i=0;i<dom.nseg;i++) {
				if(site>dom.getSegment(i).to) {
					dom1.getSegment(dom1.nseg).to=dom.getSegment(i).to;
					dom1.getSegment(dom1.nseg).from=dom.getSegment(i).from;
					dom1.nseg++;
					dom1.size+=(dom.getSegment(i).to - dom.getSegment(i).from + 1);
				}
				else if(site<dom.getSegment(i).from) {
					dom2.getSegment(dom2.nseg).to=dom.getSegment(i).to;
					dom2.getSegment(dom2.nseg).from=dom.getSegment(i).from;
					dom2.nseg++;
					dom2.size+=(dom.getSegment(i).to - dom.getSegment(i).from + 1);
				}
				else if(site>dom.getSegment(i).from &&
						site<dom.getSegment(i).to) {
					dom1.getSegment(dom1.nseg).from=dom.getSegment(i).from;
					dom1.getSegment(dom1.nseg).to=site-1;
					dom1.nseg++;
					dom1.size+=(site-dom.getSegment(i).from);
					dom2.getSegment(dom2.nseg).to=dom.getSegment(i).to;
					dom2.getSegment(dom2.nseg).from=site;
					dom2.nseg++;
					dom2.size+=(dom.getSegment(i).to-site+1);
				}
			}
		}
		else if(val.site2>0) { /* double cut */
			for(i=0;i<dom.nseg;i++) {
				if(site>dom.getSegment(i).to||val.site2<dom.getSegment(i).from) {
					dom1.getSegment(dom1.nseg).to=dom.getSegment(i).to;
					dom1.getSegment(dom1.nseg).from=dom.getSegment(i).from;
					dom1.nseg++;
					dom1.size+=(dom.getSegment(i).to - dom.getSegment(i).from + 1);
				}
				else if(site<dom.getSegment(i).from&&val.site2>dom.getSegment(i).to) {
					dom2.getSegment(dom1.nseg).to=dom.getSegment(i).to;
					dom2.getSegment(dom1.nseg).from=dom.getSegment(i).from;
					dom2.nseg++;
					dom2.size+=(dom.getSegment(i).to - dom.getSegment(i).from + 1);
				}
				else if(site>dom.getSegment(i).from &&
						site<dom.getSegment(i).to) {
					dom1.getSegment(dom1.nseg).to=site;
					dom1.getSegment(dom1.nseg).from=dom.getSegment(i).from;
					dom1.size+=(dom1.getSegment(dom1.nseg).to - dom1.getSegment(dom1.nseg).from + 1);
					dom1.nseg++;
					dom2.getSegment(dom2.nseg).from=site+1;
					if(val.site2>dom.getSegment(i).from &&
							val.site2<dom.getSegment(i).to) {
						dom2.getSegment(dom2.nseg).to=val.site2-1;
						dom2.size+=(dom2.getSegment(dom2.nseg).to - dom2.getSegment(dom2.nseg).from + 1);
						dom2.nseg++;
						dom1.getSegment(dom1.nseg).from = val.site2;
						dom1.getSegment(dom1.nseg).to = dom.getSegment(i).to;
						dom1.size+=(dom1.getSegment(dom1.nseg).to - dom1.getSegment(dom1.nseg).from + 1);
						dom1.nseg++;
					}
					else {
						dom2.getSegment(dom2.nseg).to=dom.getSegment(i).to;
						dom2.size+=(dom2.getSegment(dom2.nseg).to - dom2.getSegment(dom2.nseg).from + 1);
						dom2.nseg++;
					}
				}
				else if(val.site2>dom.getSegment(i).from &&
						val.site2<dom.getSegment(i).to) {
					dom2.getSegment(dom2.nseg).to=val.site2-1;
					dom2.getSegment(dom2.nseg).from=dom.getSegment(i).from;
					dom2.size+=(dom2.getSegment(dom2.nseg).to - dom2.getSegment(dom2.nseg).from + 1);
					dom2.nseg++;
					dom1.getSegment(dom1.nseg).from=val.site2;
					dom1.getSegment(dom1.nseg).to = dom.getSegment(i).to;
					dom1.size+=(dom1.getSegment(dom1.nseg).to - dom1.getSegment(dom1.nseg).from + 1);
					dom1.nseg++;
				}
			}
		}
		if(verbose) 
			System.out.println(String.format("  CUTR dom1 ...  nseg %d",dom1.nseg));
		
		if ( verbose)
		for(i=0;i<dom1.nseg;i++) 
			System.out.println(String.format("	F ... from %d to %d",dom1.getSegment(i).from,dom1.getSegment(i).to));
		
		cutDomain(dom1, cut_sites, pdpMatrix);
		
		if(verbose) 
			System.out.println(String.format("  C ... cutr dom2: nseg %d",dom2.nseg));
		if(verbose)
			for(i=0;i<dom2.nseg;i++) 
			 System.out.println(String.format("	F ... from %d to %d",dom2.getSegment(i).from,dom2.getSegment(i).to));
		
		cutDomain(dom2, cut_sites, pdpMatrix);
		
		//System.out.println("end of cutDomain 0 " +dom);
		//System.out.println("end of cutDomain 1 " +dom1);
		//System.out.println("end of cutDomain 2 " +dom2);
		
	}


	public List<Domain> getDomains() {
		
		return domains;
	}


	
}
