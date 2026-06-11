import sys
import sqlite3
import concurrent.futures

# Function to write a subset of table entries to a text file
def write_entries_to_file(query_pattern, start_id, end_id, file_name, db_path):
    # Connect to the SQLite3 database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Select entries in the specified range
    cursor.execute(query_pattern, (start_id, end_id))

    # Write entries to the specified file
    with open(file_name, "w") as file:
        for row in cursor.fetchall():
            file.write(str(row) + "\n")

    # Close the database connection
    conn.close()


# Function to divide tasks and run them in parallel
def distribute_entries_to_files(db_path, table_name, query_pattern, output_files):
    num_files = len(output_files)

    # Connect to the SQLite3 database to calculate splits
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Get the total number of entries in the table
    cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
    total_entries = cursor.fetchone()[0]

    # Calculate the number of entries per file (split size)
    entries_per_file = total_entries // num_files
    extra_entries = total_entries % num_files

    # Use ThreadPoolExecutor to write to files in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_files) as executor:
        futures = []
        for i, file_name in enumerate(output_files):
            start_id = i * entries_per_file + 1
            end_id = start_id + entries_per_file - 1
            if i == num_files - 1:  # Add any extra entries to the last file
                end_id += extra_entries
            futures.append(
                executor.submit(
                    write_entries_to_file,
                    query_pattern,
                    start_id,
                    end_id,
                    file_name,
                    db_path,
                )
            )

        # Wait for all tasks to complete
        concurrent.futures.wait(futures)

    # Close the database connection
    conn.close()


# Example usage
if __name__ == "__main__":
    db_path = sys.argv[1]
    table_name = sys.argv[2]
    output_files = sys.argv[3:]
    query_pattern = f"SELECT kmer FROM count_ WHERE rowid BETWEEN ? AND ?"
    distribute_entries_to_files(db_path, table_name, query_pattern, output_files)
