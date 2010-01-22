// $Id: SimpleMutableMsa.java,v 1.5 2009/01/13 19:49:30 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.List;

public class SimpleMutableMsa implements Msa {

    final List<Character>[] _data;
    final String[]          _names;

    public SimpleMutableMsa( final String[] names, final String[] seqs ) {
        _names = names;
        _data = new List[ seqs.length ];
        int row = 0;
        for( final String seq : seqs ) {
            final List<Character> l = new ArrayList<Character>( seq.length() );
            for( int i = 0; i < seq.length(); i++ ) {
                l.add( seq.charAt( i ) );
            }
            _data[ row++ ] = l;
        }
    }

    public void deleteResidueAt( final int col, final int row ) {
        //TODO
    }

    public double getGapScore() {
        //TODO
        return 0;
    }

    @Override
    public int getLength() {
        return _data[ 0 ].size();
    }

    @Override
    public int getNumberOfSequences() {
        return _data.length;
    }

    @Override
    public char getResidueAt( final int col, final int row ) {
        return _data[ row ].get( col );
    }

    public double getScore() {
        //TODO
        return 0;
    }

    @Override
    public Sequence getSequence( final int i ) {
        throw new AssertionError( "method not implemented" );
    }

    public void insertResidueAt( final int col, final int row ) {
        //TODO
    }
}
