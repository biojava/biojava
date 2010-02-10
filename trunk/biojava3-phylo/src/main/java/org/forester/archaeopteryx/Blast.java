// $Id: Blast.java,v 1.3 2009/11/03 00:01:01 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class Blast {

    public Blast() {
    }

    /**
     * Method for access REST
     * @param query
     * service name¡¢method name and parameter for exection rest
     * @return
     * execution result
     * @throws IOException
     */
    private String getResult( final String query ) throws IOException {
        final String baseURL = "http://xml.nig.ac.jp/rest/Invoke";
        final URL url = new URL( baseURL );
        final URLConnection urlc = url.openConnection();
        urlc.setDoOutput( true );
        urlc.setAllowUserInteraction( false );
        final PrintStream ps = new PrintStream( urlc.getOutputStream() );
        ps.print( query );
        ps.close();
        final BufferedReader br = new BufferedReader( new InputStreamReader( urlc.getInputStream() ) );
        final StringBuffer sb = new StringBuffer();
        String line = null;
        while ( ( line = br.readLine() ) != null ) {
            sb.append( line + "\n" );
        }
        br.close();
        return sb.toString();
    }

    void go( final String geneName ) {
        // Retrieve accession number list which has specified gene name from searchByXMLPath of ARSA. Please click here for details of ARSA.
        /*target: Sequence length is bettween 300bp and 1000bp.
        Feature key is CDS.
        Gene qualifire is same as specified gene name.*/
        String queryPath = "/ENTRY/DDBJ/division=='HUM' AND (/ENTRY/DDBJ/length>=300 AND "
                + "/ENTRY/DDBJ/length<=1000) ";
        queryPath += "AND (/ENTRY/DDBJ/feature-table/feature{/f_key = 'CDS' AND ";
        queryPath += "/f_quals/qualifier{/q_name = 'gene' AND /q_value=='" + geneName + "'}})";
        String query = "service=ARSA&method=searchByXMLPath&queryPath=" + queryPath
                + "&returnPath=/ENTRY/DDBJ/primary-accession&offset=1&count=100";
        //Execute ARSA
        String arsaResult = null;
        try {
            arsaResult = getResult( query );
        }
        catch ( final IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        final String[] arsaResultLines = arsaResult.split( "\n" );
        //Get hit count
        final int arsaResultNum = Integer.parseInt( arsaResultLines[ 0 ].replaceAll( "hitscount       =", "" ).trim() );
        //If there is no hit, print a message and exit
        if ( arsaResultNum == 0 ) {
            System.out.println( "There is no entry for gene:" + geneName );
            return;
        }
        //Retrieve DNA sequence of top hit entry by using getFASTA_DDBJEntry of GetEntry.
        //Retrieve DNA sequence of first fit.
        final String repAccession = arsaResultLines[ 2 ];
        query = "service=GetEntry&method=getFASTA_DDBJEntry&accession=" + repAccession;
        String dnaSeq = null;
        try {
            dnaSeq = getResult( query );
        }
        catch ( final IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        System.out.println( "Retrieved DNA sequence is: " + dnaSeq );
        //Execute blastn by using searchParam of Blast with step2's sequence. Specified option is -e 0.0001 -m 8 -b 50 -v 50. It means "Extract top 50 hit which E-value is more than 0.0001.". The reference databases are specified as follows. ddbjpri(primates) ddbjrod(rodents) ddbjmam(mammals) ddbjvrt(vertebrates ) ddbjinv(invertebrates).
        //Execute blastn with step3's sequence
        query = "service=Blast&method=searchParam&program=blastn&database=ddbjpri ddbjrod ddbjmam ddbjvrt "
                + "ddbjinv&query=" + dnaSeq + "&param=-m 8 -b 50 -v 50 -e 0.0001";
        String blastResult = null;
        try {
            blastResult = getResult( query );
        }
        catch ( final IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        // Extract both accession number and similarity score from BLAST result.
        // This step does not use Web API and extract the part of result or edit the result. Please click here to see the details of each column in the BLAST tab delimited format which is generated by -m 8 option.
        final String blastResultLines[] = blastResult.split( "\n" );
        final Vector<String[]> parsedBlastResult = new Vector<String[]>();
        for( final String blastResultLine : blastResultLines ) {
            final String cols[] = blastResultLine.split( "\t" );
            final String accession = cols[ 1 ].substring( 0, cols[ 1 ].indexOf( "|" ) );
            final String[] result = { accession, cols[ 2 ] };
            parsedBlastResult.add( result );
        }
        // Retrieve species name by using searchByXMLPath of ARSA. If the plural subjects whose species
        // name are same are in the result, the highest similarity score is used.
        //Retrieve species from accession number.
        final Hashtable<String, String> organismAccession = new Hashtable<String, String>();
        for( int i = 0; i < parsedBlastResult.size(); i++ ) {
            final String[] parsed = parsedBlastResult.elementAt( i );
            query = "service=ARSA&method=searchByXMLPath&queryPath=/ENTRY/DDBJ/primary-accession=='" + parsed[ 0 ]
                    + "'&returnPath=/ENTRY/DDBJ/organism&offset=1&count=100";
            String organism = null;
            try {
                organism = getResult( query );
            }
            catch ( final IOException e ) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            final String[] organismLines = organism.split( "\n" );
            organism = organismLines[ 2 ];
            //If same organism name hits, use first hit.
            if ( !organismAccession.containsKey( organism ) ) {
                organismAccession.put( organism, parsed[ 0 ] + "\t" + parsed[ 1 ] );
            }
        }
        // Print result.
        // Print Result
        System.out.println( "DDBJ entries: " + arsaResultNum );
        System.out.println( "Representative accession: " + repAccession );
        System.out.println( "Organism name\tDDBJ accession number\tSequence similarity" );
        final String[] keys = new String[ organismAccession.size() ];
        final Enumeration<String> enu = organismAccession.keys();
        int count = 0;
        while ( enu.hasMoreElements() ) {
            keys[ count ] = enu.nextElement();
            ++count;
        }
        Arrays.sort( keys );
        for( final String key : keys ) {
            System.out.println( key + "\t" + organismAccession.get( key ) );
        }
    }
}
