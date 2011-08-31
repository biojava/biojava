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
				if(site>dom.getSegmentAtPos(i).getTo()) {
					dom1.getSegmentAtPos(dom1.nseg).setTo(dom.getSegmentAtPos(i).getTo());
					dom1.getSegmentAtPos(dom1.nseg).setFrom(dom.getSegmentAtPos(i).getFrom());
					dom1.nseg++;
					dom1.size+=(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
				}
				else if(site<dom.getSegmentAtPos(i).getFrom()) {
					dom2.getSegmentAtPos(dom2.nseg).setTo(dom.getSegmentAtPos(i).getTo());
					dom2.getSegmentAtPos(dom2.nseg).setFrom(dom.getSegmentAtPos(i).getFrom());
					dom2.nseg++;
					dom2.size+=(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
				}
				else if(site>dom.getSegmentAtPos(i).getFrom() &&
						site<dom.getSegmentAtPos(i).getTo()) {
					dom1.getSegmentAtPos(dom1.nseg).setFrom(dom.getSegmentAtPos(i).getFrom());
					dom1.getSegmentAtPos(dom1.nseg).setTo(site-1);
					dom1.nseg++;
					dom1.size+=(site-dom.getSegmentAtPos(i).getFrom());
					dom2.getSegmentAtPos(dom2.nseg).setTo(dom.getSegmentAtPos(i).getTo());
					dom2.getSegmentAtPos(dom2.nseg).setFrom(site);
					dom2.nseg++;
					dom2.size+=(dom.getSegmentAtPos(i).getTo()-site+1);
				}
			}
		}
		else if(val.site2>0) { /* double cut */
			for(i=0;i<dom.nseg;i++) {
				if(site>dom.getSegmentAtPos(i).getTo()||val.site2<dom.getSegmentAtPos(i).getFrom()) {
					dom1.getSegmentAtPos(dom1.nseg).setTo(dom.getSegmentAtPos(i).getTo());
					dom1.getSegmentAtPos(dom1.nseg).setFrom(dom.getSegmentAtPos(i).getFrom());
					dom1.nseg++;
					dom1.size+=(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
				}
				else if(site<dom.getSegmentAtPos(i).getFrom()&&val.site2>dom.getSegmentAtPos(i).getTo()) {
					dom2.getSegmentAtPos(dom1.nseg).setTo(dom.getSegmentAtPos(i).getTo());
					dom2.getSegmentAtPos(dom1.nseg).setFrom(dom.getSegmentAtPos(i).getFrom());
					dom2.nseg++;
					dom2.size+=(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
				}
				else if(site>dom.getSegmentAtPos(i).getFrom() &&
						site<dom.getSegmentAtPos(i).getTo()) {
					dom1.getSegmentAtPos(dom1.nseg).setTo(site);
					dom1.getSegmentAtPos(dom1.nseg).setFrom(dom.getSegmentAtPos(i).getFrom());
					dom1.size+=(dom1.getSegmentAtPos(dom1.nseg).getTo() - dom1.getSegmentAtPos(dom1.nseg).getFrom() + 1);
					dom1.nseg++;
					dom2.getSegmentAtPos(dom2.nseg).setFrom(site+1);
					if(val.site2>dom.getSegmentAtPos(i).getFrom() &&
							val.site2<dom.getSegmentAtPos(i).getTo()) {
						dom2.getSegmentAtPos(dom2.nseg).setTo(val.site2-1);
						dom2.size+=(dom2.getSegmentAtPos(dom2.nseg).getTo() - dom2.getSegmentAtPos(dom2.nseg).getFrom() + 1);
						dom2.nseg++;
						dom1.getSegmentAtPos(dom1.nseg).setFrom( val.site2);
						dom1.getSegmentAtPos(dom1.nseg).setTo( dom.getSegmentAtPos(i).getTo());
						dom1.size+=(dom1.getSegmentAtPos(dom1.nseg).getTo() - dom1.getSegmentAtPos(dom1.nseg).getFrom() + 1);
						dom1.nseg++;
					}
					else {
						dom2.getSegmentAtPos(dom2.nseg).setTo(dom.getSegmentAtPos(i).getTo());
						dom2.size+=(dom2.getSegmentAtPos(dom2.nseg).getTo() - dom2.getSegmentAtPos(dom2.nseg).getFrom() + 1);
						dom2.nseg++;
					}
				}
				else if(val.site2>dom.getSegmentAtPos(i).getFrom() &&
						val.site2<dom.getSegmentAtPos(i).getTo()) {
					dom2.getSegmentAtPos(dom2.nseg).setTo(val.site2-1);
					dom2.getSegmentAtPos(dom2.nseg).setFrom(dom.getSegmentAtPos(i).getFrom());
					dom2.size+=(dom2.getSegmentAtPos(dom2.nseg).getTo() - dom2.getSegmentAtPos(dom2.nseg).getFrom() + 1);
					dom2.nseg++;
					dom1.getSegmentAtPos(dom1.nseg).setFrom(val.site2);
					dom1.getSegmentAtPos(dom1.nseg).setTo( dom.getSegmentAtPos(i).getTo());
					dom1.size+=(dom1.getSegmentAtPos(dom1.nseg).getTo() - dom1.getSegmentAtPos(dom1.nseg).getFrom() + 1);
					dom1.nseg++;
				}
			}
		}
		if(verbose) 
			System.out.println(String.format("  CUTR dom1 ...  nseg %d",dom1.nseg));
		
		if ( verbose)
		for(i=0;i<dom1.nseg;i++) 
			System.out.println(String.format("	F ... from %d to %d",dom1.getSegmentAtPos(i).getFrom(),dom1.getSegmentAtPos(i).getTo()));
		
		cutDomain(dom1, cut_sites, pdpMatrix);
		
		if(verbose) 
			System.out.println(String.format("  C ... cutr dom2: nseg %d",dom2.nseg));
		if(verbose)
			for(i=0;i<dom2.nseg;i++) 
			 System.out.println(String.format("	F ... from %d to %d",dom2.getSegmentAtPos(i).getFrom(),dom2.getSegmentAtPos(i).getTo()));
		
		cutDomain(dom2, cut_sites, pdpMatrix);
		
		//System.out.println("end of cutDomain 0 " +dom);
		//System.out.println("end of cutDomain 1 " +dom1);
		//System.out.println("end of cutDomain 2 " +dom2);
		
	}


	public List<Domain> getDomains() {
		
		return domains;
	}


	
}
