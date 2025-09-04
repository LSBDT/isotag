#!/usr/bin/env python3
"""
RefGet SQ.XXXX Identifier Calculator for IsoTag v1.0

Calculates RefGet sequence identifiers from FASTA files following
the official RefGet specification and GA4GH standards.

References:
- RefGet Specification: http://samtools.github.io/hts-specs/refget.html
- VRS Implementation: https://github.com/ga4gh/vrs-python
"""

import json
import re
from pathlib import Path
from typing import Dict, Set, Optional, Iterator, Tuple
import click
from vrs_compat import sha512t24u


class RefGetCalculator:
    """Calculate and manage RefGet SQ.XXXX identifiers from FASTA files"""
    
    @staticmethod
    def normalize_sequence(sequence: str) -> str:
        """Normalize sequence according to RefGet specification
        
        Args:
            sequence: Raw DNA sequence string
            
        Returns:
            Normalized sequence (uppercase, no whitespace)
        """
        # Remove all whitespace and convert to uppercase
        normalized = ''.join(sequence.upper().split())
        
        # Validate that sequence contains only valid DNA bases
        valid_chars = set('ACGTURYKMSWBDHVN')  # Include ambiguous bases
        invalid_chars = set(normalized) - valid_chars
        if invalid_chars:
            raise ValueError(f"Invalid characters in sequence: {invalid_chars}")
            
        return normalized
    
    @staticmethod 
    def calculate_refget_id(sequence: str) -> str:
        """Calculate RefGet SQ.XXXX identifier from sequence
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            RefGet identifier (e.g., SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2)
        """
        # Normalize sequence according to RefGet spec
        normalized = RefGetCalculator.normalize_sequence(sequence)
        
        # Calculate sha512t24u digest
        digest = sha512t24u(normalized.encode('ascii'))
        
        # Return with SQ. prefix
        return f"SQ.{digest}"
    
    @classmethod
    def parse_fasta_headers(cls, fasta_file: str) -> Set[str]:
        """Extract all chromosome names from FASTA headers
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            Set of chromosome names found in headers
        """
        chromosomes = set()
        
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Extract chromosome name from header
                    header = line.strip()[1:]  # Remove '>'
                    
                    # Handle different header formats
                    # >chr1 description
                    # >1 description  
                    # >chr1
                    # >NC_000001.11 Homo sapiens chromosome 1
                    
                    chrom_name = header.split()[0]  # Take first part
                    chromosomes.add(chrom_name)
        
        return chromosomes
    
    @classmethod
    def parse_fasta_sequences(cls, fasta_file: str) -> Iterator[Tuple[str, str]]:
        """Parse FASTA file and yield (chromosome, sequence) pairs
        
        Args:
            fasta_file: Path to FASTA file
            
        Yields:
            Tuple of (chromosome_name, sequence)
        """
        current_chrom = None
        current_seq = []
        
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                if line.startswith('>'):
                    # Yield previous sequence if exists
                    if current_chrom and current_seq:
                        yield current_chrom, ''.join(current_seq)
                    
                    # Start new chromosome
                    header = line[1:]  # Remove '>'
                    current_chrom = header.split()[0]  # Take first part
                    current_seq = []
                    
                elif current_chrom:  # Skip empty lines before first header
                    current_seq.append(line)
            
            # Yield last sequence
            if current_chrom and current_seq:
                yield current_chrom, ''.join(current_seq)
    
    @classmethod
    def calculate_from_fasta(cls, fasta_file: str, progress: bool = True) -> Dict[str, str]:
        """Calculate RefGet mapping for all chromosomes in FASTA
        
        Args:
            fasta_file: Path to FASTA file
            progress: Show progress messages
            
        Returns:
            Dictionary mapping chromosome names to RefGet IDs
        """
        refget_mapping = {}
        
        if progress:
            chromosomes = cls.parse_fasta_headers(fasta_file)
            click.echo(f"🧬 Processing {len(chromosomes)} chromosomes from {Path(fasta_file).name}")
        
        processed = 0
        for chrom_name, sequence in cls.parse_fasta_sequences(fasta_file):
            if progress:
                processed += 1
                click.echo(f"⏳ Processing {chrom_name} ({processed})...")
            
            # Calculate RefGet ID
            refget_id = cls.calculate_refget_id(sequence)
            refget_mapping[chrom_name] = refget_id
            
            if progress:
                click.echo(f"✅ {chrom_name} -> {refget_id}")
        
        return refget_mapping
    
    @classmethod
    def generate_all_variants(cls, base_mapping: Dict[str, str]) -> Dict[str, str]:
        """Generate all chromosome name variants for lookup
        
        Args:
            base_mapping: Basic chromosome -> RefGet mapping
            
        Returns:
            Extended mapping with all chromosome name variants
        """
        extended_mapping = {}
        
        for chrom_name, refget_id in base_mapping.items():
            # Add original name
            extended_mapping[chrom_name] = refget_id
            
            # Generate common variants
            variants = cls._generate_chromosome_variants(chrom_name)
            for variant in variants:
                extended_mapping[variant] = refget_id
        
        return extended_mapping
    
    @staticmethod
    def _generate_chromosome_variants(chrom_name: str) -> Set[str]:
        """Generate common chromosome name variants
        
        Args:
            chrom_name: Base chromosome name
            
        Returns:
            Set of chromosome name variants
        """
        variants = set()
        
        # Handle chr1 -> 1, Chr1 -> 1, CHR1 -> 1
        if chrom_name.lower().startswith('chr'):
            base_num = chrom_name[3:]
            variants.add(base_num)  # chr1 -> 1
            variants.add(f"chr{base_num}")  # Ensure lowercase
            variants.add(f"Chr{base_num}")  # Title case
            variants.add(f"CHR{base_num}")  # Uppercase
        
        # Handle 1 -> chr1, Chr1, CHR1
        elif chrom_name.isalnum():  # Simple chromosome number/letter
            variants.add(f"chr{chrom_name}")
            variants.add(f"Chr{chrom_name}")
            variants.add(f"CHR{chrom_name}")
        
        # Handle RefSeq accessions (NC_000001.11 -> 1, chr1, etc.)
        nc_match = re.match(r'NC_(\d+)\.', chrom_name)
        if nc_match:
            chrom_num = str(int(nc_match.group(1)))  # Remove leading zeros
            if chrom_num == '23':
                chrom_num = 'X'
            elif chrom_num == '24':
                chrom_num = 'Y'
            variants.add(chrom_num)
            variants.add(f"chr{chrom_num}")
            variants.add(f"Chr{chrom_num}")
            variants.add(f"CHR{chrom_num}")
        
        return variants
    
    @classmethod
    def save_refget_mapping(cls, mapping: Dict[str, str], output_file: str, 
                           genome_name: str = "unknown"):
        """Save RefGet mapping to JSON file
        
        Args:
            mapping: RefGet mapping dictionary
            output_file: Output JSON file path
            genome_name: Genome assembly name (e.g., hg38, mm39)
        """
        from datetime import datetime
        
        # Generate extended mapping with all variants
        extended_mapping = cls.generate_all_variants(mapping)
        
        data = {
            "metadata": {
                "genome": genome_name,
                "generated": datetime.now().isoformat(),
                "total_sequences": len(mapping),
                "total_mappings": len(extended_mapping),
                "description": "RefGet SQ.XXXX identifiers for chromosome sequences"
            },
            "refget_mapping": extended_mapping
        }
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2, sort_keys=True)
        
        click.echo(f"💾 Saved RefGet mapping: {output_file}")
        click.echo(f"📊 {len(mapping)} chromosomes -> {len(extended_mapping)} total mappings")
    
    @classmethod
    def load_refget_mapping(cls, json_file: str) -> Dict[str, str]:
        """Load RefGet mapping from JSON file
        
        Args:
            json_file: JSON file path
            
        Returns:
            RefGet mapping dictionary
        """
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        return data.get("refget_mapping", {})


@click.command()
@click.option('--fasta', '-f', help='Input FASTA file')
@click.option('--output', '-o', help='Output RefGet JSON mapping')
@click.option('--genome', '-g', default="unknown", help='Genome assembly name (e.g., hg38, mm39)')
@click.option('--test', is_flag=True, help='Run test calculations')
def main(fasta, output, genome, test):
    """RefGet SQ.XXXX Calculator
    
    Calculate RefGet sequence identifiers from FASTA files for use with IsoTag v1.0.
    
    Examples:
        # Calculate RefGet IDs for human genome
        python3 refget.py -f hg38.fa -o hg38-refget.json -g hg38
        
        # Calculate RefGet IDs for mouse genome  
        python3 refget.py -f mm39.fa -o mm39-refget.json -g mm39
    """
    if test:
        # Test with known examples
        test_seq = "ACGT"
        expected = "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
        result = RefGetCalculator.calculate_refget_id(test_seq)
        
        if result == expected:
            click.echo(f"✅ Test passed: {test_seq} -> {result}")
        else:
            click.echo(f"❌ Test failed: {test_seq} -> {result} (expected {expected})")
        return
    
    # Validate required parameters for non-test mode
    if not fasta or not output:
        click.echo("❌ --fasta and --output are required (use --test for testing)")
        return
    
    # Validate input file
    if not Path(fasta).exists():
        click.echo(f"❌ FASTA file not found: {fasta}")
        return
    
    click.echo(f"🧬 RefGet Calculator - Processing {fasta}")
    click.echo(f"📝 Genome: {genome}")
    
    # Calculate RefGet mapping
    refget_mapping = RefGetCalculator.calculate_from_fasta(fasta)
    
    # Save to JSON
    RefGetCalculator.save_refget_mapping(refget_mapping, output, genome)
    
    click.echo(f"\n🎯 RefGet calculation complete!")
    click.echo(f"📁 Use with IsoTag: --refget {output}")


if __name__ == "__main__":
    main()
