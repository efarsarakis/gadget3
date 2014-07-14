#!/usr/bin/perl

open(FILE, "cpu.txt");
open(OUTFILE, ">cpu_idl.txt");

$count = 0;
$compression = 50;

while($line=<FILE>)
{
    $count++;

    if($count == $compression)
    {
	chop $line;
	($first, $second)= split(":", $line);
	($time, $second)= split(",", $second);
	print OUTFILE "$time ";
    }

    for($i=0; $i<35; $i++)
    {
        $line=<FILE>;

	if($count == $compression)
	{
	    @fields = split ' ' , $line;
	    print OUTFILE "$fields[1] ";
	}
    }

    $line=<FILE>;   # empty line

    if($count == $compression)
    {
	print OUTFILE "\n";
	$count = 0;
    }
	
}


