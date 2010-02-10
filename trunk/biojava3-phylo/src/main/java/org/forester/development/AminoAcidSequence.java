// $Id: AminoAcidSequence.java,v 1.16 2009/01/13 19:49:30 cmzmasek Exp $
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

package org.forester.development;

import java.io.IOException;
import java.io.Writer;

import org.forester.phylogeny.data.PhylogenyData;

public class AminoAcidSequence implements PhylogenyData {

    private final byte[] _sequence;
    private final String _name;

    public AminoAcidSequence( final int length ) {
        _sequence = new byte[ length ];
        _name = "";
    }

    public AminoAcidSequence( final String name, final byte[] sequence ) {
        _sequence = new byte[ sequence.length ];
        _name = new String( name );
        for( int i = 0; i < sequence.length; ++i ) {
            setStateAt( i, sequence[ i ] );
        }
    }

    public AminoAcidSequence( final String name, final String sequence ) {
        _sequence = new byte[ sequence.length() ];
        _name = new String( name );
        for( int i = 0; i < sequence.length(); ++i ) {
            setResidueAt( i, sequence.charAt( i ) );
        }
    }

    public StringBuffer asSimpleText() {
        return new StringBuffer( getName() );
    }

    public StringBuffer asText() {
        return new StringBuffer( getName() ).append( ": " ).append( getSequenceAsString() );
    }

    public AminoAcidSequence copy() {
        return ( new AminoAcidSequence( getName(), _sequence ) );
    }

    public int getLength() {
        return _sequence.length;
    }

    public String getName() {
        return _name;
    }

    public char getResidueAt( final int position ) {
        return AminoAcid.getResidue( getStateAt( position ) );
    }

    public String getSequenceAsString() {
        final StringBuffer sb = new StringBuffer( getLength() );
        for( int i = 0; i < getLength(); ++i ) {
            sb.append( getResidueAt( i ) );
        }
        return sb.toString();
    }

    public byte getStateAt( final int position ) {
        return _sequence[ position ];
    }

    public boolean isEqual( final PhylogenyData data ) {
        // TODO Auto-generated method stub
        return false;
    }

    public void setResidueAt( final int position, final char residue ) {
        setStateAt( position, AminoAcid.getState( residue ) );
    }

    public void setStateAt( final int position, byte state ) {
        if ( AminoAcid.isUnassignable( state ) ) {
            state = AminoAcid.UNKNOWN_CODE;
        }
        _sequence[ position ] = state;
    }

    public StringBuffer toNHX() {
        // TODO Auto-generated method stub
        return null;
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( "<sequence type=\"aminoacid\">\n" );
        writer.write( "<name>" );
        writer.write( getName() );
        writer.write( "</name>\n" );
        writer.write( "<sequence>" );
        writer.write( getSequenceAsString() );
        writer.write( "</sequence>\n" );
        writer.write( "</sequence>\n" );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
