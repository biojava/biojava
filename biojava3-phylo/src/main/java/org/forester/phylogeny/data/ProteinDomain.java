// $Id: ProteinDomain.java,v 1.36 2008/09/24 16:42:47 cmzmasek Exp $
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class ProteinDomain implements PhylogenyData {

    final public static double CONFIDENCE_DEFAULT = 0.0;
    final public static String IDENTIFIER_DEFAULT = "";
    final private String       _name;
    final private int          _from;
    final private int          _to;
    final private String       _id;
    final private double       _confidence;

    public ProteinDomain( final String name, final int from, final int to ) {
        this( name, from, to, ProteinDomain.IDENTIFIER_DEFAULT, ProteinDomain.CONFIDENCE_DEFAULT );
    }

    public ProteinDomain( final String name, final int from, final int to, final double confidence ) {
        this( name, from, to, ProteinDomain.IDENTIFIER_DEFAULT, confidence );
    }

    public ProteinDomain( final String name, final int from, final int to, final String id ) {
        this( name, from, to, id, ProteinDomain.CONFIDENCE_DEFAULT );
    }

    public ProteinDomain( final String name, final int from, final int to, final String id, final double confidence ) {
        if ( ( from >= to ) || ( to < 0 ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain from " + from + " to " + to );
        }
        _name = name;
        _from = from;
        _to = to;
        _id = id;
        _confidence = confidence;
    }

    public StringBuffer asSimpleText() {
        return new StringBuffer( getName() );
    }

    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer( getName() );
        sb.append( " [" );
        sb.append( getLength() );
        if ( !ForesterUtil.isEmpty( getId() ) ) {
            sb.append( " " );
            sb.append( getId() );
        }
        if ( getConfidence() != CONFIDENCE_DEFAULT ) {
            sb.append( " " );
            sb.append( getConfidence() );
        }
        sb.append( "]" );
        return sb;
    }

    public PhylogenyData copy() {
        if ( getId() == null ) {
            return new ProteinDomain( new String( getName() ), getFrom(), getTo(), getConfidence() );
        }
        return new ProteinDomain( new String( getName() ), getFrom(), getTo(), new String( getId() ), getConfidence() );
    }

    public double getConfidence() {
        return _confidence;
    }

    public int getFrom() {
        return _from;
    }

    public String getId() {
        return _id;
    }

    public int getLength() {
        return ( getTo() - getFrom() + 1 );
    }

    public String getName() {
        return _name;
    }

    public int getTo() {
        return _to;
    }

    public boolean isEqual( final PhylogenyData protein_domain ) {
        if ( protein_domain == null ) {
            return false;
        }
        if ( !( protein_domain instanceof ProteinDomain ) ) {
            return false;
        }
        else if ( ( ( ProteinDomain ) protein_domain ).getLength() != getLength() ) {
            return false;
        }
        else if ( !( ( ProteinDomain ) protein_domain ).getName().equals( getName() ) ) {
            return false;
        }
        return true;
    }

    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        if ( getId() != null ) {
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_DOMAIN,
                                          PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_FROM,
                                          getFrom() + "",
                                          PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_TO,
                                          getTo() + "",
                                          PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_CONFIDENCE,
                                          getConfidence() + "",
                                          PhyloXmlMapping.IDENTIFIER,
                                          getId() );
        }
        else {
            PhylogenyDataUtil.appendOpen( writer,
                                          PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_DOMAIN,
                                          PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_FROM,
                                          getFrom() + "",
                                          PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_TO,
                                          getTo() + "",
                                          PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_CONFIDENCE,
                                          getConfidence() + "" );
        }
        writer.write( getName() );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_DOMAIN );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
