// $Id: DirectedBinaryDomainCombination.java,v 1.4 2008/03/09 00:11:13 cmzmasek
// Exp $
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

package org.forester.surfacing;

public class AdjactantDirectedBinaryDomainCombination extends BasicBinaryDomainCombination {

    public AdjactantDirectedBinaryDomainCombination( final DomainId n_terminal, final DomainId c_terminal ) {
        super();
        if ( ( n_terminal == null ) || ( c_terminal == null ) ) {
            throw new IllegalArgumentException( "attempt to create binary domain combination using null" );
        }
        _id_0 = n_terminal;
        _id_1 = c_terminal;
    }

    public AdjactantDirectedBinaryDomainCombination( final String n_terminal, final String c_terminal ) {
        this( new DomainId( n_terminal ), new DomainId( c_terminal ) );
    }

    public static AdjactantDirectedBinaryDomainCombination createInstance( final String ids ) {
        if ( ids.indexOf( BinaryDomainCombination.SEPARATOR ) < 1 ) {
            throw new IllegalArgumentException( "Unexpected format for binary domain combination [" + ids + "]" );
        }
        final String[] ids_ary = ids.split( BinaryDomainCombination.SEPARATOR );
        if ( ids_ary.length != 2 ) {
            throw new IllegalArgumentException( "Unexpected format for binary domain combination [" + ids + "]" );
        }
        return new AdjactantDirectedBinaryDomainCombination( ids_ary[ 0 ], ids_ary[ 1 ] );
    }
}
