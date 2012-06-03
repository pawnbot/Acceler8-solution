#include <malloc.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <sched.h>
#include <string.h>

#include "processor_op.h"

#define CPU_INFO "/proc/cpuinfo"

void get_proc_info (int* cpu_count, int* numa_nodes, cpu_set_t* affinity)
{
  int iter;
  char buff[2048];
  cpu_set_t active;
  FILE* info_input;
  char const* pos;
  char const* str_proc = "processor";
  char const* str_node = "physical id";
  int cpu_id, node_id;
  
  CPU_ZERO (&active);
  sched_getaffinity(0, sizeof(active), &active);
  *cpu_count = 0;
  
  for (iter = 0; iter < CPU_SETSIZE; ++iter)
  {
    if (CPU_ISSET (iter, &active))
    {
      ++(*cpu_count);
    }
  }
  
  return;                    
  info_input = fopen(CPU_INFO, "r");
  if (info_input != 0 && affinity != 0) {
    *numa_nodes = 1;
    node_id = -1;
    cpu_id = -1;
    while (fgets(buff, 2048, info_input)) {
      if (strncmp(buff, str_proc, sizeof(str_proc)) == 0)
      {
        pos = strchr(buff + sizeof(str_proc), ':');
        if (pos)
          cpu_id = atoi(pos + 1);
      }
      if (strncmp(buff, str_node, sizeof(str_node)) == 0)
      {
        pos = strchr(buff + sizeof(str_node), ':');
        if (pos)
          node_id = atoi(pos + 1);
      }
      if (node_id + 1 > *numa_nodes)
        *numa_nodes = node_id + 1;
        
      if (node_id >= 0 && cpu_id >= 0) {
        CPU_SET(cpu_id, affinity + node_id);
        cpu_id = -1; node_id = -1;
      } 
    }
  }
  fclose(info_input);
}
